#!/bin/bash

# Set working directory and input paths
cwd="/Users/xd14188/Desktop/UoB/1kg_hg38/"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
dbSNPv157_pos="/Users/xd14188/Desktop/UoB/tools/genomeRef/dbSNPv157_hg38/dbsnp157_biallelic_rs_hg38_UCSC.tab"
populations=("EUR")

vcf_pop_dir="${cwd}/vcfPops"
# Create output directories if missing
mkdir -p "${cwd}/tmp_vcf_sampleName"
mkdir -p "${cwd}/notInVcfSample"
mkdir -p "${cwd}/bedFiles_maf001"
mkdir -p "${vcf_pop_dir}"





# seprate by population and filter with MAF > 0.01
for chr in {1..22,x}; do
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
        vcf_pop_out="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf.gz"
        vcf_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList.tab"
        flt_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList_filtered.tab"
        vcf_flt_out="${vcf_pop_dir}/chr${chr}_${pop}_only_filtered_SNP.vc.gz"
        bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

        # Check if sample list exists
        if [[ ! -f "$sample_list" ]]; then
            echo "❌ Sample list not found: $sample_list"
            continue
        fi

        # Identify missing samples and log them
        comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

        # Extract samples and filter by dbSNP
        echo "Extracting and filtering for chr$chr $pop..."
        bcftools view -S "$sample_list" --force-samples "$vcf_file" -Ou | \
        bcftools +fill-tags -Ou -- -t MAF | \
        bcftools view -i 'MAF>0.01' -Oz -o "$vcf_pop_out"
        bcftools index "$vcf_pop_out"
        echo "--- Done: chr$chr $pop  population Sep and MAF Cal ---"

        # make the "$vcf_pop_out" sample variants list
        bcftools query -f '%CHROM\t%POS\n' "$vcf_pop_out" > "$vcf_var_out"
        # compare the variants in the vcf file with the dbSNP v157 hg38 position tab file
        LC_ALL=C awk 'NR==FNR{a[$1 FS $2];next} ($1 FS $2) in a' \
        "$vcf_var_out" "$dbSNPv157_pos" > "$flt_var_out"
        echo "--- Done: chr$chr $pop  SNPs check in dbSNPv157 ---"
        
        # filter "$vcf_out" variants with the dbSNP v157 hg38 matched variants list "$flt_var_out"
        bcftools view -R "$flt_var_out" -Oz -o "$vcf_flt_out" "$vcf_pop_out"
        echo "--- Done: chr$chr $pop ---"
        
        # Clean up intermediate VCF
        # rm -f "$vcf_unzipped"

        echo "-------------------------"
    done
done