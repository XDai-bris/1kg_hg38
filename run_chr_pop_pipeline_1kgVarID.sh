#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
# Set working directory and paths
cwd="/user/home/xd14188/repo/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
dbsnp_dir="/user/home/xd14188/repo/data/genomeRef/cpra_rsID" 
r_script="${cwd}/filter_multiallelic_by_freq.R"

populations=("EUR" "AFR" "EAS" "SAS" "AMR")

# chroms=({1..22} X)
chroms=({1..22})

# Output base directories
mkdir -p "${cwd}/logs" "${cwd}/bedFiles_maf001"

# === FUNCTION ===
process_chr_pop() {
    chr="$1"
    pop="$2"

    echo "🔄 chr${chr} | ${pop}"

    log_file="${cwd}/logs/chr${chr}_${pop}.log"
    vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    keep_file="${sample_dir}/${pop}.fam"
    dbsnp_file="${dbsnp_dir}/chr${chr}_dbsnp.tsv"

    tmpdir="${cwd}/tmp/chr${chr}_${pop}"
    outdir="${cwd}/bedFiles_maf001/chr${chr}_${pop}"
    unmatched_dir="${cwd}/filterVar"
    multiallelic_dir="${cwd}/fltMultiallelics"
    unmapped_rsID_dir="${cwd}/fltUnMatchedRsID"

    mkdir -p "$tmpdir" "$outdir" "$unmatched_dir" "$multiallelic_dir" "$unmapped_rsID_dir"

    out_prefix="$tmpdir/out"
    cleaned_prefix="$tmpdir/out_cleaned"
    renamed_bim="$tmpdir/renamed.bim"
    unmatched_rsID_bim="$unmapped_rsID_dir/chr${chr}_${pop}_unmatched_rsID.bim"
    multiallelic_bim="$multiallelic_dir/chr${chr}_${pop}_multiallelic.bim"
    final_prefix="${outdir}/final_output"

    {
        echo "=== Processing chr${chr} | ${pop} ==="

        # --- Step 0: Input Validation ---
        [[ ! -f "$vcf_file" ]] && echo "❌ Missing VCF: $vcf_file" && exit 1
        [[ ! -f "$keep_file" ]] && echo "❌ Missing sample list: $keep_file" && exit 1
        [[ ! -f "$dbsnp_file" ]] && echo "❌ dbSNP file: $dbsnp_file" && exit 1

        # --- Step 1: Convert VCF to PLINK ---
        echo "[STEP 1] Convert VCF to PLINK with MAF ≥ 0.01..."
        plink2 --vcf "$vcf_file" --keep "$keep_file" --maf 0.01 \
            --make-bed --out "$out_prefix" --silent

        # --- Step 2: Compute Allele Frequencies ---
        echo "[STEP 2] Compute allele frequencies..."
        afreq_file="$tmpdir/afreq_chr${chr}_${pop}.afreq"
        plink2 --bfile "$out_prefix" --freq --out "${afreq_file%.afreq}" --silent

        # --- Step 3: Filter multiallelic variants via R ---
        echo "[STEP 3] Run R script to filter multiallelic variants..."
        filtered_prefix="$tmpdir/filtered_output_chr${chr}_${pop}"
        Rscript "$r_script" "${out_prefix}.bim" "$afreq_file" "$filtered_prefix"

        # --- Step 4: Keep only highest-frequency variants ---
        echo "[STEP 4] Filter .bed/.bim/.fam using keep list..."
        plink2 --bfile "$out_prefix" \
            --extract "${filtered_prefix}_keep.txt" \
            --make-bed --out "$cleaned_prefix" --silent

        # --- Step 5: Save removed multiallelics (BIM only) ---
        echo "[STEP 5] Save removed multiallelic variants (BIM only)..."
        awk 'NR==FNR{a[$1]; next} ($2 in a)' \
            "${filtered_prefix}_remove.txt" "${out_prefix}.bim" \
            > "$multiallelic_bim"

        # --- Step 6: Rename variants in BIM using dbSNP mapping ---
        echo "[STEP 6] Rename variants in BIM..."
        awk '
            BEGIN { OFS="\t" }
            FNR == NR { map[$1] = $2; next }
            {
                key = $1":"$4":"$6":"$5
                if (key in map) {
                    $2 = map[key]
                    print > "'"$renamed_bim"'"
                } else {
                    print > "'"$unmatched_rsID_bim"'"
                }
            }
        ' "$dbsnp_file" "${cleaned_prefix}.bim"

        # --- Step 7: Generate Final Output ---
        echo "[STEP 7] Generating final PLINK file with rsID..."
        comm -23 <(cut -f2 "${cleaned_prefix}.bim" | sort) \
                 <(cut -f2 "$unmatched_rsID_bim" | sort) > "$tmpdir/keep_final.txt"

        plink2 --bfile "$cleaned_prefix" --extract "$tmpdir/keep_final.txt" \
            --make-bed --out "$final_prefix" --silent

        cp "$renamed_bim" "${final_prefix}.bim"

        # --- Step 8: Verify Final Output with --freq ---
        echo "[STEP 8] Verifying final PLINK file with --freq..."
        plink2 --bfile "$final_prefix" --freq \
            --out "${final_prefix}_freq" --silent

        # --- Final Summary ---
        unmapped_count=$(wc -l < "$unmatched_rsID_bim")
        multiallelic_count=$(wc -l < "$multiallelic_bim")

        echo ""
        echo "✅ Done: chr${chr} ${pop}"
        echo "→ Final BED/BIM/FAM with rsIDs: $final_prefix.*"
        echo "→ Removed multiallelics (BIM only): $multiallelic_bim (n = $multiallelic_count)"
        echo "→ Unmatched variants (no rsID match): $unmatched_rsID_bim (n = $unmapped_count)"
        echo "→ Intermediate files in: $tmpdir/"
    } &> "$log_file"
}

export -f process_chr_pop
export cwd vcf_dir sample_dir dbsnp_dir plink2 r_script

# === RUN PARALLEL ===
echo "🚀 Launching parallel jobs..."
parallel -j 12 process_chr_pop ::: "${chroms[@]}" ::: "${populations[@]}"
echo "✅ All jobs submitted."
