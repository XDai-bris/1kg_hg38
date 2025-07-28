# ----MAKE BED FROM MERGED VCF FILE ---- #
#!/bin/bash

set -euo pipefail

# Config
base="/user/home/xd14188/repo/1kg_hg38"
vcf_dir="${base}/mergedPopVcf"
out_dir="${base}/mergedPopBed"
mkdir -p "$out_dir"

# List of populations
populations=("EUR" "AFR" "EAS" "SAS" "AMR")

# Convert function
convert_to_bed() {
  pop="$1"
  input_vcf="${vcf_dir}/${pop}_merged.vcf.gz"
  output_prefix="${out_dir}/${pop}_plink"

  if [[ ! -f "$input_vcf" ]]; then
    echo "❌ Missing input VCF: $input_vcf"
    return 1
  fi

  echo "🔄 Converting $input_vcf to PLINK binary format..."

  plink2 --vcf "$input_vcf" \
         --make-bed \
         --out "$output_prefix" \
         --threads 1 \
         --silent

  echo "✅ Done: ${output_prefix}.bed"
}
export -f convert_to_bed
export vcf_dir out_dir

# Run in parallel (adjust -j for concurrency)
parallel -j 5 convert_to_bed ::: "${populations[@]}"
