#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
cwd="/user/home/xd14188/repo/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
dbsnp_file="/user/home/xd14188/repo/data/genomeRef/cpra_rsID/chrX_dbsnp.tsv"
r_script="${cwd}/filter_multiallelic_by_freq.R"

populations=("EUR" "AFR" "EAS" "SAS" "AMR")
chr="X"

# === MAIN ===
for pop in "${populations[@]}"; do
    echo "🔄 chrX | ${pop}"

    log_file="${cwd}/logs/chrX_${pop}.log"
    vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    keep_file="${sample_dir}/${pop}.fam"

    tmpdir="${cwd}/tmp/chrX_${pop}"
    outdir="${cwd}/bedFiles_maf001/chrX_${pop}"
    unmatched_dir="${cwd}/filterVar"
    multiallelic_dir="${cwd}/fltMultiallelics"
    unmapped_rsID_dir="${cwd}/fltUnMatchedRsID"

    mkdir -p "$tmpdir" "$outdir" "$unmatched_dir" "$multiallelic_dir" "$unmapped_rsID_dir"

    out_prefix="$tmpdir/out"
    cleaned_prefix="$tmpdir/out_cleaned"
    renamed_bim="$tmpdir/renamed.bim"
    unmatched_rsID_bim="$unmapped_rsID_dir/chrX_${pop}_unmatched_rsID.bim"
    multiallelic_bim="$multiallelic_dir/chrX_${pop}_multiallelic.bim"
    final_prefix="${outdir}/final_output"

    {
        echo "=== Processing chrX | ${pop} ==="

        [[ ! -f "$vcf_file" ]] && echo "❌ Missing VCF: $vcf_file" && exit 1
        [[ ! -f "$keep_file" ]] && echo "❌ Missing sample list: $keep_file" && exit 1
        [[ ! -f "$dbsnp_file" ]] && echo "❌ dbSNP file: $dbsnp_file" && exit 1

        # --- Step 1: Convert VCF to PLINK ---
        echo "[STEP 1] Convert VCF to PLINK with MAF ≥ 0.01..."
        plink2 --vcf "$vcf_file" --keep "$keep_file" --maf 0.01 \
            --make-bed --out "$out_prefix" --silent

        # --- Step 1.5: Impute missing variant IDs as CHR:POS:REF:ALT ---
        echo "[STEP 1.5] Imputing missing variant IDs as CHR:POS:REF:ALT..."
        awk 'BEGIN{OFS="\t"} { $2 = $1 ":" $4 ":" $6 ":" $5; print }' "${out_prefix}.bim" > "${out_prefix}.bim.tmp"
        mv "${out_prefix}.bim.tmp" "${out_prefix}.bim"

        echo "[STEP 2] Compute allele frequencies..."
        afreq_file="$tmpdir/afreq_chrX_${pop}.afreq"
        plink2 --bfile "$out_prefix" --freq --out "${afreq_file%.afreq}" --silent

        echo "[STEP 3] Filter multiallelic variants via R..."
        filtered_prefix="$tmpdir/filtered_output_chrX_${pop}"
        Rscript "$r_script" "${out_prefix}.bim" "$afreq_file" "$filtered_prefix"

        echo "[STEP 4] Filter .bed/.bim/.fam using keep list..."
        plink2 --bfile "$out_prefix" \
            --extract "${filtered_prefix}_keep.txt" \
            --make-bed --out "$cleaned_prefix" --silent

        echo "[STEP 5] Save removed multiallelic variants (BIM only)..."
        awk 'NR==FNR{a[$1]; next} ($2 in a)' \
            "${filtered_prefix}_remove.txt" "${out_prefix}.bim" \
            > "$multiallelic_bim"

        echo "[STEP 6] Rename variants in BIM (chrX-aware)..."
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

        echo "[STEP 7] Generate final PLINK file with rsID..."
        comm -23 <(cut -f2 "${cleaned_prefix}.bim" | sort) \
                 <(cut -f2 "$unmatched_rsID_bim" | sort) > "$tmpdir/keep_final.txt"

        plink2 --bfile "$cleaned_prefix" --extract "$tmpdir/keep_final.txt" \
            --make-bed --out "$final_prefix" --silent

        cp "$renamed_bim" "${final_prefix}.bim"

        echo "[STEP 8] Verify final PLINK with --freq"
        plink2 --bfile "$final_prefix" --freq \
            --out "${final_prefix}_freq" --silent

        unmapped_count=$(wc -l < "$unmatched_rsID_bim")
        multiallelic_count=$(wc -l < "$multiallelic_bim")

        echo ""
        echo "✅ Done: chrX ${pop}"
        echo "→ Final: $final_prefix.*"
        echo "→ Unmatched rsIDs: $unmatched_rsID_bim (n = $unmapped_count)"
        echo "→ Multiallelic removed: $multiallelic_bim (n = $multiallelic_count)"
        echo "→ Intermediate: $tmpdir/"
    } &> "$log_file"
done
