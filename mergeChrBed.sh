#!/bin/bash
# Merge per-chromosome BED sets into genome-wide sets per population on SLURM.
# - One srun task per population (no GNU Parallel).
# - Per-task logs written to OUT_DIR.
# - Uses plink2 --threads=${SLURM_CPUS_PER_TASK}.
# - Handles PLINK2 a4.3: use --pmerge-list with 3 columns for BED/BIM/FAM.

set -euo pipefail
export PS1="${PS1-}"   # avoid PS1 unbound in profile scripts under set -u

# ------------------ CONFIG ------------------
BASE_DIR="/user/home/xd14188/repo/1kg_hg38/bedFiles_maf001"
OUT_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome"
MERGE_SUMMARY="${OUT_DIR}/variants_pop.txt"

populations=("AFR" "AMR" "EAS" "EUR" "SAS")
chroms=({1..22} X)

T="${SLURM_CPUS_PER_TASK:-8}"
[[ "${DEBUG:-0}" == "1" ]] && set -x

# ------------------ HELPERS ------------------
load_modules() {
  set +u
  if command -v module >/dev/null 2>&1; then
    source /etc/profile 2>/dev/null || true
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    module load plink2/2.00a4.3-netlib-lapack-qbk6
  fi
  set -u
  command -v plink2 >/dev/null 2>&1 || { echo "❌ plink2 not found on PATH" >&2; exit 127; }
}

merge_one_pop () {
  local pop="$1"
  echo "🔄 Processing population: ${pop}"
  export PS1="${PS1-}"
  load_modules

  export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMEXPR_NUM_THREADS=$T

  mkdir -p "${OUT_DIR}"
  local prefixes_file="${OUT_DIR}/${pop}_prefixes.txt"
  : > "${prefixes_file}"

  # Build list of per-chr prefixes we actually have
  local missing=0
  for chr in "${chroms[@]}"; do
    local pref="${BASE_DIR}/chr${chr}_${pop}/final_output"
    if [[ -f "${pref}.bed" && -f "${pref}.bim" && -f "${pref}.fam" ]]; then
      echo "${pref}" >> "${prefixes_file}"
    else
      echo "⚠️  Missing BED/BIM/FAM for ${pref}" >&2
      ((missing++)) || true
    fi
  done

  local nsets
  nsets=$(grep -vc '^[[:space:]]*$' "${prefixes_file}" || echo 0)
  if (( nsets == 0 )); then
    echo "❌ ${pop}: no per-chrom files found; skipping." >&2
    return 1
  fi

  local out="${OUT_DIR}/${pop}_genome"
  echo "📦 Merging ${nsets} chromosome sets for ${pop} (threads=${T})..."

  if (( nsets == 1 )); then
    # Single fileset → just rewrite to genome-wide prefix
    local only; only=$(head -n1 "${prefixes_file}")
    plink2 --bfile "${only}" --make-bed --out "${out}" --silent --threads "${T}"
  else
    # PLINK2 a4.3 requires 3-column pmerge list for BED inputs
    local pmerge3="${OUT_DIR}/${pop}_pmerge_list_3col.txt"
    : > "${pmerge3}"
    while read -r pref; do
      echo "${pref}.bed ${pref}.bim ${pref}.fam" >> "${pmerge3}"
    done < "${prefixes_file}"

    # Merge ALL sets in one go (no base --bfile); output stays BED-mode
    plink2 --pmerge-list "${pmerge3}" \
           --make-bed --out "${out}" --silent --threads "${T}"
  fi

  local nvar; nvar=$(wc -l < "${out}.bim")
  echo -e "${pop}\t${nvar}" >> "${MERGE_SUMMARY}"

  echo "📊 Allele frequencies for ${pop}..."
  plink2 --bfile "${out}" --freq --out "${out}_freq" --silent --threads "${T}"

  echo "✅ ${pop} done. Variants: ${nvar}. Missing-chr sets during scan: ${missing}."
}

# ------------------ ENTRYPOINT ------------------
mkdir -p "${OUT_DIR}"

# If called with a population arg, run just that pop
if [[ $# -eq 1 ]]; then
  merge_one_pop "$1"
  exit 0
fi

# Launcher mode
: > "${MERGE_SUMMARY}"
echo "🚀 Launching merge jobs with srun (threads per task = ${T})..."

set +u  # don't leak nounset into subshells that source profiles
pids=()
for p in "${populations[@]}"; do
  srun --exclusive -N1 -n1 --cpus-per-task="${T}" --mpi=none \
       --output="${OUT_DIR}/merge_${p}-%x-%j-%N-%t-%r.out" \
       --error="${OUT_DIR}/merge_${p}-%x-%j-%N-%t-%r.err" \
       env -u SHELLOPTS bash -c "/user/home/xd14188/repo/1kg_hg38/mergeChrBed.sh ${p}" &
  pids+=($!)
done

fail=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then fail=1; fi
done
set -u

if (( fail )); then
  echo "❌ One or more population merges failed. Check ${OUT_DIR}/merge_*.err" >&2
  exit 1
fi

echo "📋 Variant counts saved to: ${MERGE_SUMMARY}"
echo "✅ All populations merged and allele frequencies computed."