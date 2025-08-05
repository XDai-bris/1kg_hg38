#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
PLINK2="plink2"
PLINK="plink"  # Assuming PLINK 1.9 is in your PATH
BASE_DIR="/user/home/xd14188/repo/1kg_hg38/bedFiles_maf001"
OUT_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome"
MERGE_SUMMARY="$OUT_DIR/variants_pop.txt"

# Set population(s) to process
populations=("AFR")
chroms=$(seq 1 22)

mkdir -p "$OUT_DIR"
> "$MERGE_SUMMARY"  # Clear old summary

for pop in "${populations[@]}"; do
    echo "🔄 Processing population: $pop"

    pmerge_file="$OUT_DIR/${pop}_pmerge_list.txt"
    > "$pmerge_file"

    # === Step 1: Create merge list ===
    for chr in $chroms; do
        prefix="${BASE_DIR}/chr${chr}_${pop}/final_output"
        if [[ -f "${prefix}.bed" ]]; then
            echo "$prefix" >> "$pmerge_file"
        else
            echo "⚠️  Missing: ${prefix}.bed"
        fi
    done

    output_prefix="${OUT_DIR}/${pop}_genome"
    merge_list="$pmerge_file"
    current_suffix=""

    # === Step 2: Merge with auto-clean loop ===
    echo "📦 Merging chromosomes for $pop (auto-clean)..."

    $PLINK --merge-list "$merge_list" --make-bed --out "${output_prefix}${current_suffix}" --allow-no-sex --silent || true

    while [[ -f "${output_prefix}${current_suffix}-merge.missnp" ]]; do
        echo "❌ Merge failed — cleaning problematic variants..."

        cleaned_merge_list="${merge_list%.txt}_cleaned${current_suffix}.txt"
        > "$cleaned_merge_list"

        while read -r bfile; do
            cleaned_bfile="${bfile}_cleaned${current_suffix}"
            $PLINK --bfile "$bfile" \
                   --exclude "${output_prefix}${current_suffix}-merge.missnp" \
                   --make-bed \
                   --out "$cleaned_bfile" \
                   --silent
            echo "$cleaned_bfile" >> "$cleaned_merge_list"
        done < "$merge_list"

        current_suffix="${current_suffix}_clean"
        output_prefix_clean="${OUT_DIR}/${pop}_genome${current_suffix}"
        merge_list="$cleaned_merge_list"

        $PLINK --merge-list "$merge_list" --make-bed --out "$output_prefix_clean" --allow-no-sex --silent || true
    done

    final_output="${OUT_DIR}/${pop}_genome${current_suffix}"
    variant_count=$(wc -l < "${final_output}.bim")
    echo -e "${pop}\t${variant_count}" >> "$MERGE_SUMMARY"
    echo "✅ Done: $pop (Variants: $variant_count)"

    # === Step 3: Generate allele frequencies ===
    echo "📊 Generating allele frequencies for $pop..."
    $PLINK2 --bfile "$final_output" \
            --freq \
            --out "${final_output}_freq" \
            --silent
done

echo ""
echo "📋 Variant counts saved to: $MERGE_SUMMARY"
echo "✅ All populations merged and allele frequencies computed."
