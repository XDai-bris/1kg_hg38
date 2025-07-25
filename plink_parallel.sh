#!/bin/bash

set -euo pipefail

# Config
cwd="/Users/xd14188/Desktop/UoB/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
dbSNPv157_pos="/Users/xd14188/Desktop/UoB/tools/genomeRef/dbSNPv157_hg38/dbsnp157_biallelic_rs_hg38_UCSC.tab"
populations=("EUR" "AFR" "EAS" "SAS" "AMR")
vcf_pop_dir="${cwd}/vcfPops"

# Create directories
mkdir -p "${cwd}/tmp_vcf_sampleName" "${cwd}/notInVcfSample" "${cwd}/bedFiles_maf001" "${vcf_pop_dir}"

# ---- FUNCTION ---- #
process_chr_pop() {
    chr="$1"
    pop="$2"

    echo "=== Processing chr${chr} | ${pop} ==="

    vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    vcf_sample_list="${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"

    sample_list="${sample_dir}/sampleName_${pop}.txt"
    not_in_vcf="${cwd}/notInVcfSample/chr${chr}_notInVcfSample_${pop}.txt"
    vcf_pop_out="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf.gz"
    vcf_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList.tab"
    flt_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList_filtered.tab"
    vcf_flt_out="${vcf_pop_dir}/chr${chr}_${pop}_only_filtered_SNP.vc.gz"
    bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

    # Check input files
    [[ ! -f "$vcf_file" ]] && echo "❌ Missing VCF: $vcf_file" && return
    [[ ! -f "$sample_list" ]] && echo "❌ Missing sample list: $sample_list" && return

    # Extract VCF sample names (if not already done)
    if [[ ! -f "$vcf_sample_list" ]]; then
        bcftools query -l "$vcf_file" > "$vcf_sample_list"
    fi

    # Track missing samples
    comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

    # Filter samples, compute MAF, and extract variants with MAF > 0.01
    echo "-> Filtering samples and computing MAF..."
    bcftools view -S "$sample_list" --force-samples -m2 -M2 -v snps "$vcf_file" -Ou | \
    bcftools +fill-tags -Ou -- -t MAF | \
    bcftools view -i 'MAF>0.01' -Oz -o "$vcf_pop_out"

    bcftools index "$vcf_pop_out"

    # Count number of variants after filtering
    variant_count=$(bcftools view -H "$vcf_pop_out" | wc -l)
    echo "🧬 Variants after filtering: $variant_count"

    # Create list of variants
    bcftools query -f '%CHROM\t%POS\n' "$vcf_pop_out" > "$vcf_var_out"

    # Intersect with dbSNP positions
    LC_ALL=C awk 'NR==FNR{a[$1 FS $2]; next} ($1 FS $2) in a' "$vcf_var_out" "$dbSNPv157_pos" > "$flt_var_out"

    # Filter final VCF with dbSNP intersected variants
    bcftools view -R "$flt_var_out" -Oz -o "$vcf_flt_out" "$vcf_pop_out"

    bcftools index "$vcf_flt_out"

    echo "✅ Done: chr${chr} ${pop}"
}
export -f process_chr_pop
export cwd vcf_dir sample_dir dbSNPv157_pos vcf_pop_dir

# ---- PARALLEL EXECUTION ---- #
echo "🚀 Starting parallel job..."

# Define chromosomes (1–22 and X)
chroms=$(seq 1 22) && chroms+=" X"

# Launch jobs in parallel (adjust -j for concurrency level)
parallel -j 2 process_chr_pop ::: $chroms ::: "${populations[@]}"

echo "🎉 All done!"
