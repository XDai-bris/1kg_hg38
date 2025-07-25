#!/bin/bash

# Set working directory
cwd="/Users/xd14188/Desktop/UoB/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
populations=("EUR" "AFR" "EAS" "SAS" "AMR")
# here to choose to use dbSNP v157 hg38 or dbSNP v137 hg38 with biallelic and rs ID valid SNPs
dbSNPv157_hg38_dir="/Users/xd14188/Desktop/UoB/tools/genomeRef/dbSNPv157_hg38/dbsnp157_biallelic_rs_hg38.bed"

# Create output directories
mkdir -p "${cwd}/tmp_vcf_sampleName"
mkdir -p "${cwd}/notInVcfSample"
mkdir -p "${cwd}/bedFiles_maf001"
mkdir -p "${cwd}/vcfPops"
# Extract sample names from each VCF ahead of time
for chr in {1..22}; do
  vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  bcftools query -l "$vcf_file" > "${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"
done

# Initialize summary table
summary_file="${cwd}/tmp_vcf_sampleName/sample_name_summary.tsv"
echo -e "Chromosome\tMatchWithChr1" > "$summary_file"

# Use chr1 as reference
ref_file="${cwd}/tmp_vcf_sampleName/chr1_vcf_samples.txt"

# Compare other chromosomes to chr1
for chr in {1..22}; do
  target_file="${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"
  
  if cmp -s "$ref_file" "$target_file"; then
    echo -e "chr${chr}\tYES" >> "$summary_file"
  else
    echo -e "chr${chr}\tNO" >> "$summary_file"
  fi
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
  vcf_out="${cwd}/vcfPops/chr${chr}_${pop}_only.vcf.gz"
  vcf_unzipped="${cwd}/vcfPops/chr${chr}_${pop}_only.vcf"
  bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

  if [[ ! -f "$sample_list" || ! -f "$vcf_file" ]]; then
    echo "❌ Missing file for chr${chr} $pop ❌"
    return
  fi

  comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"
    # Extract population samples (ignore missing with --force-samples)
    # Use bcftools to filter VCF by sample list and bed file from dbSNP v137 hg38
    # Note: Adjust the path to your dbSNP BED file as needed
  echo "Extracting samples for chr$chr $pop..."
  bcftools view -S "$sample_list" --force-samples "$vcf_file" \
  | bcftools view -R "$dbSNPv157_hg38_dir" \
  -Oz -o "$vcf_out"
  gunzip -k "$vcf_out"
  plink --vcf "$vcf_unzipped" --make-bed --maf 0.01 --out "$bed_prefix"
  rm "$vcf_unzipped"  # Clean up unzipped VCF
}

export -f process_chr_pop

# Run the jobs in parallel: 4 concurrent tasks
parallel -j 3 process_chr_pop ::: {1..1} ::: EUR AFR EAS SAS AMR
