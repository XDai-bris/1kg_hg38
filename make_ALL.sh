#!/usr/bin/env bash
# make_ALL.sh — PLINK 1.9 union merge of AFR/AMR/EAS/EUR/SAS, then PGEN + .afreq
set -euo pipefail
export PS1="${PS1-}"

# -------- PATHS --------
IN_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome"
OUT_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome_ALL"
TMP_DIR="${OUT_DIR}/_tmp"
POPS=(AFR AMR EAS EUR SAS)   # AFR is base/template

# Threads
T="${SLURM_CPUS_PER_TASK:-8}"
[[ "${DEBUG:-0}" == "1" ]] && set -x

# -------- MODULES --------
have(){ command -v "$1" >/dev/null 2>&1; }
load_modules() {
  set +u
  if command -v module >/dev/null 2>&1; then
    source /etc/profile 2>/dev/null || true
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    module load plink/1.9-beta6.27-openblas-ibxp
    module load plink2/2.00a4.3-netlib-lapack-qbk6
  fi
  set -u
  have plink  || { echo "❌ plink 1.9 not found"; exit 127; }
  have plink2 || { echo "❌ plink2 not found"; exit 127; }
}
load_modules

export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMEXPR_NUM_THREADS=$T
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# Helpers
dedup_prefix(){ echo "${TMP_DIR}/$1_dedup"; }
harm_prefix(){  echo "${TMP_DIR}/$1_std";   }

# -------- SANITY CHECK INPUTS --------
for P in "${POPS[@]}"; do
  for ext in bed bim fam; do
    [[ -f "${IN_DIR}/${P}_genome.${ext}" ]] || { echo "❌ Missing ${IN_DIR}/${P}_genome.${ext}" >&2; exit 2; }
  done
done

# -------- STEP 0: Bullet-proof per-pop de-dup by BIM (IDs) --------
echo "🧼 De-duplicating variant IDs from BIM per population (keep first occurrence)..."
for P in "${POPS[@]}"; do
  # find duplicate IDs purely from BIM
  cut -f2 "${IN_DIR}/${P}_genome.bim" | sort | uniq -d > "${TMP_DIR}/${P}_dups.txt" || true
  N=$(wc -l < "${TMP_DIR}/${P}_dups.txt" || echo 0)
  if (( N > 0 )); then
    echo "  • ${P}: removing ${N} duplicate IDs"
    plink \
      --bfile "${IN_DIR}/${P}_genome" \
      --exclude "${TMP_DIR}/${P}_dups.txt" \
      --make-bed \
      --out "$(dedup_prefix ${P})" \
      >/dev/null
  else
    echo "  • ${P}: no duplicates"
    cp "${IN_DIR}/${P}_genome.bed" "$(dedup_prefix ${P}).bed"
    cp "${IN_DIR}/${P}_genome.bim" "$(dedup_prefix ${P}).bim"
    cp "${IN_DIR}/${P}_genome.fam" "$(dedup_prefix ${P}).fam"
  fi
  # verify zero dups now
  if cut -f2 "$(dedup_prefix ${P}).bim" | sort | uniq -d | grep -q .; then
    echo "❌ ${P}: duplicates still present after dedup — aborting." >&2
    exit 3
  fi
done

# -------- STEP 1: Harmonize each pop to AFR (flip if needed; drop irreconcilables) --------
echo "🧭 Harmonizing allele orientation vs AFR..."
cp "$(dedup_prefix AFR).bed" "${TMP_DIR}/AFR_std.bed"
cp "$(dedup_prefix AFR).bim" "${TMP_DIR}/AFR_std.bim"
cp "$(dedup_prefix AFR).fam" "${TMP_DIR}/AFR_std.fam"

for P in AMR EAS EUR SAS; do
  echo "  → ${P}: align to AFR"
  plink \
    --bfile  "${TMP_DIR}/AFR_std" \
    --bmerge "$(dedup_prefix ${P}).bed" "$(dedup_prefix ${P}).bim" "$(dedup_prefix ${P}).fam" \
    --make-bed \
    --out    "${TMP_DIR}/AFR_vs_${P}_try1" \
    >/dev/null 2>&1 || true

  if [[ -f "${TMP_DIR}/AFR_vs_${P}_try1-merge.missnp" ]]; then
    echo "    • flipping ${P} sites listed in missnp..."
    plink \
      --bfile "$(dedup_prefix ${P})" \
      --flip  "${TMP_DIR}/AFR_vs_${P}_try1-merge.missnp" \
      --make-bed \
      --out   "${TMP_DIR}/${P}_flipped1" \
      >/dev/null

    plink \
      --bfile  "${TMP_DIR}/AFR_std" \
      --bmerge "${TMP_DIR}/${P}_flipped1.bed" "${TMP_DIR}/${P}_flipped1.bim" "${TMP_DIR}/${P}_flipped1.fam" \
      --make-bed \
      --out    "${TMP_DIR}/AFR_vs_${P}_try2" \
      >/dev/null 2>&1 || true

    if [[ -f "${TMP_DIR}/AFR_vs_${P}_try2-merge.missnp" ]]; then
      echo "    • excluding stubborn discordant sites for ${P}..."
      plink \
        --bfile   "${TMP_DIR}/${P}_flipped1" \
        --exclude "${TMP_DIR}/AFR_vs_${P}_try2-merge.missnp" \
        --make-bed \
        --out     "$(harm_prefix ${P})" \
        >/dev/null
    else
      mv "${TMP_DIR}/${P}_flipped1".bed "$(harm_prefix ${P})".bed
      mv "${TMP_DIR}/${P}_flipped1".bim "$(harm_prefix ${P})".bim
      mv "${TMP_DIR}/${P}_flipped1".fam "$(harm_prefix ${P})".fam
    fi
  else
    cp "$(dedup_prefix ${P}).bed" "$(harm_prefix ${P})".bed
    cp "$(dedup_prefix ${P}).bim" "$(harm_prefix ${P})".bim
    cp "$(dedup_prefix ${P}).fam" "$(harm_prefix ${P})".fam
  fi
done

# -------- STEP 2: PLINK 1.9 union merge (auto-exclude discordant until clean) --------
MERGE_LIST="${TMP_DIR}/merge_list.txt"
: > "${MERGE_LIST}"
for P in AMR EAS EUR SAS; do
  echo "$(harm_prefix ${P}).bed $(harm_prefix ${P}).bim $(harm_prefix ${P}).fam" >> "${MERGE_LIST}"
done

echo "🧩 Running union merge..."
ATTEMPT=1
BASE_PREFIX="${TMP_DIR}/ALL_merge_base_${ATTEMPT}"
cp "${TMP_DIR}/AFR_std".bed "${BASE_PREFIX}.bed"
cp "${TMP_DIR}/AFR_std".bim "${BASE_PREFIX}.bim"
cp "${TMP_DIR}/AFR_std".fam "${BASE_PREFIX}.fam"

while : ; do
  plink \
    --bfile      "${BASE_PREFIX}" \
    --merge-list "${MERGE_LIST}" \
    --make-bed \
    --out        "${TMP_DIR}/ALL_merge_out_${ATTEMPT}" \
    >/dev/null 2>&1 || true

  MISS="${TMP_DIR}/ALL_merge_out_${ATTEMPT}-merge.missnp"
  if [[ -f "${MISS}" ]]; then
    echo "    • attempt ${ATTEMPT}: excluding $(wc -l < "${MISS}") discordant variants and retrying..."
    plink --bfile "${BASE_PREFIX}" --exclude "${MISS}" --make-bed --out "${TMP_DIR}/base_filtered_${ATTEMPT}" >/dev/null
    NEW_LIST="${TMP_DIR}/merge_list_attempt_${ATTEMPT}.txt"
    : > "${NEW_LIST}"
    for P in AMR EAS EUR SAS; do
      plink --bfile "$(harm_prefix ${P})" --exclude "${MISS}" --make-bed \
            --out "${TMP_DIR}/${P}_filtered_${ATTEMPT}" >/dev/null
      echo "${TMP_DIR}/${P}_filtered_${ATTEMPT}.bed ${TMP_DIR}/${P}_filtered_${ATTEMPT}.bim ${TMP_DIR}/${P}_filtered_${ATTEMPT}.fam" >> "${NEW_LIST}"
    done
    MERGE_LIST="${NEW_LIST}"
    ATTEMPT=$((ATTEMPT+1))
    BASE_PREFIX="${TMP_DIR}/base_filtered_$((ATTEMPT-1))"
    continue
  fi

  mv "${TMP_DIR}/ALL_merge_out_${ATTEMPT}.bed" "${OUT_DIR}/ALL_genome.bed"
  mv "${TMP_DIR}/ALL_merge_out_${ATTEMPT}.bim" "${OUT_DIR}/ALL_genome.bim"
  mv "${TMP_DIR}/ALL_merge_out_${ATTEMPT}.fam" "${OUT_DIR}/ALL_genome.fam"
  echo "✅ PLINK 1.9 merge completed in ${ATTEMPT} attempt(s)."
  break
done

# -------- STEP 3: PGEN + AFREQ with PLINK 2 --------
echo "🧱 Writing PGEN for ALL (plink2)..."
plink2 --bfile "${OUT_DIR}/ALL_genome" --make-pgen --out "${OUT_DIR}/ALL_genome" --threads "${T}" --silent

echo "📊 Computing allele frequencies (.afreq) (plink2)..."
plink2 --bfile "${OUT_DIR}/ALL_genome" --freq --out "${OUT_DIR}/ALL_genome_freq" --threads "${T}" --silent

echo "🎉 Done."
echo "Outputs in ${OUT_DIR}:"
echo "  • ALL_genome.bed  ALL_genome.bim  ALL_genome.fam"
echo "  • ALL_genome.pgen ALL_genome.pvar ALL_genome.psam"
echo "  • ALL_genome_freq.afreq (+ .log)"
echo "  • Working files in ${TMP_DIR} (safe to delete)"