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

# Extract sample names from each VCF ahead of time
for chr in {1..22}; do
  vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  bcftools query -l "$vcf_file" > "${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"
done

# Export the processing function
process_chr_pop() {
  chr="$1"
  pop="$2"
  cwd="/Users/xd14188/Desktop/UoB/1kg_hg38"
  vcf_dir="${cwd}/vcfFiles"
  sample_dir="${cwd}/sampleName"

  sample_list="${sample_dir}/sampleName_${pop}.txt"
  vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  vcf_sample_list="${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"
  not_in_vcf="${cwd}/notInVcfSample/chr${chr}_notInVcfSample_${pop}.txt"
  vcf_out="${vcf_dir}/chr${chr}_${pop}_only.vcf.gz"
  bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

  if [[ ! -f "$sample_list" || ! -f "$vcf_file" ]]; then
    echo "❌ Missing file for chr${chr} $pop"
    return
  fi

  comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"
  bcftools view -S "$sample_list" --force-samples -Oz -o "$vcf_out" "$vcf_file"
  gunzip -f "$vcf_out"
  vcf_unzipped="${vcf_out%.gz}"
  plink --vcf "$vcf_unzipped" --make-bed --maf 0.01 --out "$bed_prefix"
}

export -f process_chr_pop

# Run the jobs in parallel: 4 concurrent tasks
parallel -j 4 process_chr_pop ::: {1..22} ::: EUR AFR EAS SAS AMR
