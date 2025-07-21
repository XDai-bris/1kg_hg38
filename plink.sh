#!/bin/bash

# Set working directory
cwd="/Users/xd14188/Desktop/UoB/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
populations=("EUR" "AFR" "EAS" "SAS" "AMR")

# Create output directories
mkdir -p "${cwd}/tmp_vcf_sampleName"
mkdir -p "${cwd}/notInVcfSample"
mkdir -p "${cwd}/bedFiles_maf001"

# Loop through chromosomes 1 to 22
for chr in {1..2}; do
    echo "=== Processing chromosome $chr ==="
    vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    vcf_sample_list="${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"

    # Check if VCF file exists
    if [[ ! -f "$vcf_file" ]]; then
        echo "❌ VCF not found for chr$chr: $vcf_file"
        continue
    fi

    # Extract sample names from VCF
    bcftools query -l "$vcf_file" > "$vcf_sample_list"

    # Loop through each population
    for pop in "${populations[@]}"; do
        echo "--- Processing $pop for chr$chr ---"

        sample_list="${sample_dir}/sampleName_${pop}.txt"
        not_in_vcf="${cwd}/notInVcfSample/chr${chr}_notInVcfSample_${pop}.txt"
        vcf_out="${vcf_dir}/chr${chr}_${pop}_only.vcf.gz"
        bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

        # Check if sample list exists
        if [[ ! -f "$sample_list" ]]; then
            echo "❌ Sample list not found: $sample_list"
            continue
        fi

        # Identify population samples not in VCF and save list
        comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

        # Extract population samples (ignore missing with --force-samples)
        bcftools view -S "$sample_list" --force-samples -Oz -o "$vcf_out" "$vcf_file"

        # Unzip VCF for PLINK v1.9
        gunzip -f "$vcf_out"
        vcf_unzipped="${vcf_out%.gz}"

        # Convert to BED and filter by MAF > 0.01
        plink --vcf "$vcf_unzipped" --make-bed --maf 0.01 --out "$bed_prefix"

        echo "--- Done: chr$chr $pop ---"
    done
done

echo "✅ All chromosomes and populations processed."
