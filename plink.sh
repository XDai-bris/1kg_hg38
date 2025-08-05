#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
PLINK2="/Users/xd14188/Desktop/UoB/tools/plink2_mac_arm64_20250707/plink2"
VCF="/Users/xd14188/Desktop/UoB/1kg_hg38/vcfFiles/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
KEEP="/Users/xd14188/Desktop/UoB/1kg_hg38/sampleName/AFR.fam"
DBSNPFILE="/Users/xd14188/Desktop/UoB/tools/genomeRef/dbSNPv157_hg38/cpra_rsID/chr21_dbsnp.tsv"
R_SCRIPT="/Users/xd14188/Desktop/UoB/1kg_hg38/filter_multiallelic_by_freq.R"

TMPDIR="tmp"
OUTDIR="output"
UNMATCHED_DIR="filterVar"
mkdir -p "$TMPDIR" "$OUTDIR" "$UNMATCHED_DIR"

OUT_PREFIX="$TMPDIR/test_output"
CLEANED_PREFIX="$TMPDIR/test_output_cleaned"
RENAMED_PREFIX="$OUTDIR/final_output"
RENAMED_BIM="$TMPDIR/renamed.bim"
UNMATCHED_BIM="${UNMATCHED_DIR}/noMatchedRsID.bim"

# === STEP 1: Convert VCF to PLINK ===
echo "[STEP 1] Convert VCF to PLINK with MAF ≥ 0.01..."
"$PLINK2" \
  --vcf "$VCF" \
  --keep "$KEEP" \
  --maf 0.01 \
  --make-bed \
  --out "$OUT_PREFIX" \
  --silent

# === STEP 2: Compute Allele Frequencies ===
echo "[STEP 2] Compute allele frequencies..."
"$PLINK2" \
  --bfile "$OUT_PREFIX" \
  --freq \
  --out "$TMPDIR/afreq" \
  --silent

# === STEP 3: Filter multiallelic variants via R ===
echo "[STEP 3] Run R script to filter multiallelic variants..."
Rscript "$R_SCRIPT" "${OUT_PREFIX}.bim" "$TMPDIR/afreq.afreq" "$TMPDIR/filtered_output"

# === STEP 4: Keep only highest-frequency variants ===
echo "[STEP 4] Filter .bed/.bim/.fam using keep list..."
"$PLINK2" \
  --bfile "$OUT_PREFIX" \
  --extract "$TMPDIR/filtered_output_keep.txt" \
  --make-bed \
  --out "$CLEANED_PREFIX" \
  --silent

# === STEP 5: Save removed multiallelics as BIM only ===
echo "[STEP 5] Save removed multiallelic variants (BIM only)..."
awk 'NR==FNR{a[$1]; next} ($2 in a)' \
  "$TMPDIR/filtered_output_remove.txt" "${OUT_PREFIX}.bim" \
  > "$UNMATCHED_DIR/removed_multiallelic.bim"

# === STEP 6: Check dbSNP file ===
echo "[STEP 6] Checking dbSNP CHR:POS:REF:ALT → rsID file..."
if [[ ! -f "$DBSNPFILE" ]]; then
  echo "❌ dbSNP file not found: $DBSNPFILE"
  exit 1
fi

# === STEP 7: Rename Variants in BIM ===
echo "[STEP 7] Rename variants in BIM..."
BIM="${CLEANED_PREFIX}.bim"
awk '
  BEGIN { OFS="\t" }
  FNR == NR { map[$1] = $2; next }
  {
    key = $1":"$4":"$6":"$5  # CHR:POS:REF:ALT
    if (key in map) {
      $2 = map[key]
      print > "'"$RENAMED_BIM"'"
    } else {
      print > "'"$UNMATCHED_BIM"'"
    }
  }
' "$DBSNPFILE" "$BIM"

# === STEP 8: Generate Final Output ===
echo "[STEP 8] Generating final PLINK file with rsID..."

FINAL_OUT_DIR="output"
FINAL_OUT="${FINAL_OUT_DIR}/final_output"
mkdir -p "$FINAL_OUT_DIR"

# Get variant IDs to keep (i.e., those NOT in unmatched)
comm -23 <(cut -f2 "$CLEANED_PREFIX.bim" | sort) <(cut -f2 "$UNMATCHED_BIM" | sort) > "$TMPDIR/keep_final.txt"

# Filter the cleaned PLINK set
"$PLINK2" \
  --bfile "$CLEANED_PREFIX" \
  --extract "$TMPDIR/keep_final.txt" \
  --make-bed \
  --out "$FINAL_OUT" \
  --silent

# Overwrite .bim with rsID-renamed .bim
cp "$RENAMED_BIM" "$FINAL_OUT.bim"

# === STEP 9: Verify Final Output with --freq ===
echo "[STEP 9] Verifying final PLINK file with --freq..."
"$PLINK2" \
  --bfile "$FINAL_OUT" \
  --freq \
  --out "${FINAL_OUT}_freq" \
  --silent

# === DONE ===
UNMATCHED_COUNT=$(wc -l < "$UNMATCHED_BIM")
echo ""
echo "✅ Done"
echo "→ Final Output Saved:"
echo "→ Final BED/BIM/FAM with rsIDs: $FINAL_OUT.*"
echo "→ Unmatched variants: $UNMATCHED_BIM (n = $UNMATCHED_COUNT)"
echo "→ Removed multiallelics (BIM only): $UNMATCHED_DIR/removed_multiallelic.bim"
echo "→ Intermediate files saved in: $TMPDIR/"
