#!/bin/bash

set -euo pipefail

base_dir="/user/home/xd14188/repo/1kg_hg38"
vcf_dir="${base_dir}/vcfFiles"
vcf_pop_dir="${base_dir}/vcfPops"
out_dir="${base_dir}/tmp_check_varNum"
mkdir -p "$out_dir"

populations=("EUR" "AFR" "EAS" "SAS" "AMR")

for pop in "${populations[@]}"; do
  echo "🔍 Processing $pop..."

  summary_file="${out_dir}/variant_counts_${pop}_summary.txt"
  echo -e "chr\traw_lines\traw_file\tpop_lines\tpop_file\tfiltered_lines\tfiltered_file" > "$summary_file"

  for chr in {1..22} X; do
    raw_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    pop_file="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf.gz"
    flt_file="${vcf_pop_dir}/chr${chr}_${pop}_only_filtered_SNP.vcf.gz"

    raw_count="NA"
    pop_count="NA"
    flt_count="NA"

    [[ -f "$raw_file" ]] && raw_count=$(wc -l < "$raw_file")
    [[ -f "$pop_file" ]] && pop_count=$(wc -l < "$pop_file")
    [[ -f "$flt_file" ]] && flt_count=$(wc -l < "$flt_file")

    echo -e "chr${chr}\t$raw_count\t$(basename "$raw_file")\t$pop_count\t$(basename "$pop_file")\t$flt_count\t$(basename "$flt_file")" >> "$summary_file"
  done

  echo "✅ Saved: $summary_file"
done

echo "🏁 All population summaries complete in: $out_dir"



sed -i 's|/user/home/xd14188/repo/1kg_hg38/vcfFiles/1kGP_high_coverage_Illumina.| |g' "${out_dir}/"*
sed -i 's|1kGP_high_coverage_Illumina.| |g' "${out_dir}/"*
sed -i 's|.filtered.SNV_INDEL_SV_phased_panel.vcf.gz| |g' "${out_dir}/"*

echo "✅ Done. See: ${out_dir}/*.summary.txt"
