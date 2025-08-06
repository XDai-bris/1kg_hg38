#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
PLINK2="/Users/xd14188/Desktop/UoB/tools/plink2_mac_arm64_20250707/plink2"
BASE_DIR="/user/home/xd14188/repo/1kg_hg38/bedFiles_maf001"
OUT_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome"
MERGE_SUMMARY="$OUT_DIR/variants_pop.txt"

populations=("AFR" "AMR" "EAS" "EUR" "SAS")
chroms=({1..22} X)

mkdir -p "$OUT_DIR"
> "$MERGE_SUMMARY"  # Clear old summary

for pop in "${populations[@]}"; do
    echo "🔄 Processing population: $pop"

    pmerge_file="$OUT_DIR/${pop}_pmerge_list.txt"
    > "$pmerge_file"

    # === Step 1: Create pmerge list ===
    for chr in $chroms; do
        prefix="${BASE_DIR}/chr${chr}_${pop}/final_output"
        if [[ -f "${prefix}.bed" ]]; then
            echo "$prefix" >> "$pmerge_file"
        else
            echo "⚠️  Missing: ${prefix}.bed"
        fi
    done

    output_prefix="${OUT_DIR}/${pop}_genome"

    # === Step 2: Merge PLINK sets ===
    echo "📦 Merging chromosomes for $pop..."
    plink \
        --merge-list "$pmerge_file" \
        --make-bed \
        --out "$output_prefix" \
        --silent

    # === Step 3: Count variants ===
    variant_count=$(wc -l < "${output_prefix}.bim")
    echo -e "${pop}\t${variant_count}" >> "$MERGE_SUMMARY"
    echo "✅ Done: $pop (Variants: $variant_count)"

    # === Step 4: Generate allele frequencies ===
    echo "📊 Generating allele frequencies for $pop..."
    plink2 \
        --bfile "$output_prefix" \
        --freq \
        --out "${output_prefix}_freq" \
        --silent
done

echo ""
echo "📋 Variant counts saved to: $MERGE_SUMMARY"
echo "✅ All populations merged and allele frequencies computed."
