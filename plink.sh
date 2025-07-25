#!/bin/bash

# Set working directory and input paths
cwd="/Users/xd14188/Desktop/UoB/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
vcf_pop_dir="${cwd}/vcfPops"

# Choose dbSNP version (must be a BED with biallelic rsID SNVs)
dbSNPv157_bed="/Users/xd14188/Desktop/UoB/tools/genomeRef/dbSNPv157_hg38/dbsnp157_biallelic_rs_hg38.bed"

# Populations
populations=("EUR" "AFR" "EAS" "SAS" "AMR")

# Create output directories if missing
mkdir -p "${cwd}/tmp_vcf_sampleName"
mkdir -p "${cwd}/notInVcfSample"
mkdir -p "${cwd}/bedFiles_maf001"
mkdir -p "${vcf_pop_dir}"

# Loop through chromosomes 1 to 22
for chr in {1..22}; do
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
        vcf_out="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf.gz"
        vcf_unzipped="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf"
        bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

        # Check if sample list exists
        if [[ ! -f "$sample_list" ]]; then
            echo "❌ Sample list not found: $sample_list"
            continue
        fi

        # Identify missing samples and log them
        comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

        # Remove any existing conflicting output
        rm -f "$vcf_out" "$vcf_unzipped"

        # Extract samples and filter by dbSNP
        echo "Extracting and filtering for chr$chr $pop..."
        bcftools view -S "$sample_list" --force-samples "$vcf_file" \
        | bcftools view -R "$dbSNPv157_bed" -Oz -o "$vcf_out"

        # Unzip the VCF for PLINK (plink v1.9 doesn't accept .vcf.gz)
        gunzip -f "$vcf_out"

        # Convert to PLINK and filter by MAF > 0.01
        plink --vcf "$vcf_unzipped" --make-bed --maf 0.01 --out "$bed_prefix"

        # Clean up intermediate VCF
        rm -f "$vcf_unzipped"

        echo "--- Done: chr$chr $pop ---"
    done
done

echo "✅ All chromosomes and populations processed successfully."
