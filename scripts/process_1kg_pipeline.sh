#!/bin/bash
# ============================================================================
# Script: 1KGP HG38 VCF Processing Pipeline
# Purpose: Extract population-specific SNPs from 1000 Genomes VCFs,
#          merge them by population, and convert them to PLINK format.
# Author: [Xiaoyang Dai]
# Date: [28 Jul 2025]
# ============================================================================

set -euo pipefail

# --------------------------[ CONFIGURATION ]-------------------------------- #
CWD="/user/home/xd14188/repo/1kg_hg38"
VCF_DIR="${CWD}/vcfFiles"
SAMPLE_DIR="${CWD}/sampleName"
DBSNP_V157="${CWD}/../data/genomeRef/dbsnp157_biallelic_rs_hg38_UCSC.tab"
POPULATIONS=("EUR" "AFR" "EAS" "SAS" "AMR")
VCF_POP_DIR="${CWD}/vcfPops"
TMP_DIR="${CWD}/tmp_vcf_sampleName"
LOG_DIR="${VCF_POP_DIR}/logs"
NOT_IN_VCF_DIR="${CWD}/notInVcfSample"

mkdir -p "$TMP_DIR" "$NOT_IN_VCF_DIR" "$LOG_DIR"

# ----------------------[ FUNCTION: PROCESS CHR + POP ]---------------------- #
process_chr_pop() {
    local chr="$1"
    local pop="$2"
    local log_file="${LOG_DIR}/chr${chr}_${pop}.log"

    {
        echo "=== Processing chr${chr} | ${pop} ==="
        
        local vcf_file="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        local vcf_sample_list="${TMP_DIR}/chr${chr}_vcf_samples.txt"
        local sample_list="${SAMPLE_DIR}/sampleName_${pop}.txt"
        local not_in_vcf="${NOT_IN_VCF_DIR}/chr${chr}_notInVcfSample_${pop}.txt"
        local vcf_pop_out="${VCF_POP_DIR}/chr${chr}_${pop}_only.vcf.gz"
        local vcf_var_out="${VCF_POP_DIR}/chr${chr}_${pop}_only_varList.tab"
        local flt_var_out="${VCF_POP_DIR}/chr${chr}_${pop}_only_varList_filtered.tab"
        local vcf_flt_out="${VCF_POP_DIR}/chr${chr}_${pop}_only_filtered_SNP.vcf.gz"
        
        [[ ! -f "$vcf_file" ]] && echo "❌ Missing VCF: $vcf_file" && exit 1
        [[ ! -f "$sample_list" ]] && echo "❌ Missing sample list: $sample_list" && exit 1

        [[ ! -f "$vcf_sample_list" ]] && bcftools query -l "$vcf_file" > "$vcf_sample_list"

        comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

        echo "-> Filtering samples and computing MAF..."
        bcftools view -S "$sample_list" --force-samples -m2 -M2 -v snps "$vcf_file" -Ou | \
        bcftools +fill-tags -Ou -- -t MAF | \
        bcftools view -i 'MAF>0.01' -Oz -o "$vcf_pop_out"

        bcftools index "$vcf_pop_out"

        local variant_count
        variant_count=$(bcftools view -H "$vcf_pop_out" | wc -l)
        echo "🧬 Variants after filtering: $variant_count"

        bcftools query -f '%CHROM\t%POS\n' "$vcf_pop_out" > "$vcf_var_out"

        LC_ALL=C awk 'NR==FNR{a[$1 FS $2]; next} ($1 FS $2) in a' "$vcf_var_out" "$DBSNP_V157" > "$flt_var_out"

        bcftools view -R "$flt_var_out" -Oz -o "$vcf_flt_out" "$vcf_pop_out"
        bcftools index "$vcf_flt_out"

        echo "✅ Done: chr${chr} ${pop}"
    } &> "$log_file"
}
export -f process_chr_pop
export CWD VCF_DIR SAMPLE_DIR DBSNP_V157 VCF_POP_DIR

# ---------------------------[ PARALLEL RUN ]-------------------------------- #
echo "🚀 Starting parallel processing..."
CHROMS=$(echo {1..22} X)
parallel -j 8 process_chr_pop ::: $CHROMS ::: "${POPULATIONS[@]}"

# --------------------------[ MERGE PER POPULATION ]------------------------- #
OUTDIR="${CWD}/mergedPopVcf"
mkdir -p "$OUTDIR"

merge_one_pop() {
    local pop="$1"
    echo "🔄 Merging $pop..."
    local vcf_files=()
    for chr in {1..22} X; do
        local vcf="${VCF_POP_DIR}/chr${chr}_${pop}_only_filtered_SNP.vcf.gz"
        [[ -f "$vcf" ]] && vcf_files+=("$vcf")
    done
    [[ ${#vcf_files[@]} -eq 0 ]] && echo "⚠️ No files for $pop" && return

    local out="${OUTDIR}/${pop}_merged.vcf.gz"
    bcftools concat -Oz -o "$out" "${vcf_files[@]}"
    bcftools index --csi "$out"
    echo "✅ Done: $out"
}
export -f merge_one_pop
parallel -j 5 merge_one_pop ::: "${POPULATIONS[@]}" 2>&1 | tee "${OUTDIR}/merge_all.log"

# --------------------------[ CONVERT TO PLINK BED ]------------------------- #
PLINK_OUT="${CWD}/mergedPopBed"
mkdir -p "$PLINK_OUT"

convert_to_bed() {
    local pop="$1"
    local input_vcf="${OUTDIR}/${pop}_merged.vcf.gz"
    local output_prefix="${PLINK_OUT}/${pop}_plink"

    [[ ! -f "$input_vcf" ]] && echo "❌ Missing input VCF: $input_vcf" && return 1

    echo "🔄 Converting $input_vcf to PLINK binary format..."
    plink2 --vcf "$input_vcf" \
           --make-bed \
           --out "$output_prefix" \
           --threads 1 \
           --silent
    echo "✅ Done: ${output_prefix}.bed"
}
export -f convert_to_bed
export OUTDIR PLINK_OUT

parallel -j 5 convert_to_bed ::: "${POPULATIONS[@]}"
