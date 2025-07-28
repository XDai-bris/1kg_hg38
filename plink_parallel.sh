#!/bin/bash
# Set strict error handling
set -euo pipefail

# Config
cwd="/user/home/xd14188/repo/1kg_hg38"
vcf_dir="${cwd}/vcfFiles"
sample_dir="${cwd}/sampleName"
dbSNPv157_pos="/user/home/xd14188/repo/data/genomeRef/dbsnp157_biallelic_rs_hg38_UCSC.tab"
populations=("EUR" "AFR" "EAS" "SAS" "AMR")
vcf_pop_dir="${cwd}/vcfPops"

# Create necessary directories
mkdir -p "${cwd}/tmp_vcf_sampleName" "${cwd}/notInVcfSample" "${cwd}/bedFiles_maf001" "${vcf_pop_dir}/logs"

# ---- FUNCTION ---- #
process_chr_pop() {
    chr="$1"
    pop="$2"

    log_file="${vcf_pop_dir}/logs/chr${chr}_${pop}.log"
    {
        echo "=== Processing chr${chr} | ${pop} ==="

        vcf_file="${vcf_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        vcf_sample_list="${cwd}/tmp_vcf_sampleName/chr${chr}_vcf_samples.txt"

        sample_list="${sample_dir}/sampleName_${pop}.txt"
        not_in_vcf="${cwd}/notInVcfSample/chr${chr}_notInVcfSample_${pop}.txt"
        vcf_pop_out="${vcf_pop_dir}/chr${chr}_${pop}_only.vcf.gz"
        vcf_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList.tab"
        flt_var_out="${vcf_pop_dir}/chr${chr}_${pop}_only_varList_filtered.tab"
        vcf_flt_out="${vcf_pop_dir}/chr${chr}_${pop}_only_filtered_SNP.vcf.gz"
        bed_prefix="${cwd}/bedFiles_maf001/chr${chr}_${pop}_MAF01"

        [[ ! -f "$vcf_file" ]] && echo "❌ Missing VCF: $vcf_file" && exit 1
        [[ ! -f "$sample_list" ]] && echo "❌ Missing sample list: $sample_list" && exit 1

        if [[ ! -f "$vcf_sample_list" ]]; then
            bcftools query -l "$vcf_file" > "$vcf_sample_list"
        fi

        comm -23 <(sort "$sample_list") <(sort "$vcf_sample_list") > "$not_in_vcf"

        echo "-> Filtering samples and computing MAF..."
        bcftools view -S "$sample_list" --force-samples -m2 -M2 -v snps "$vcf_file" -Ou | \
        bcftools +fill-tags -Ou -- -t MAF | \
        bcftools view -i 'MAF>0.01' -Oz -o "$vcf_pop_out"

        bcftools index "$vcf_pop_out"

        variant_count=$(bcftools view -H "$vcf_pop_out" | wc -l)
        echo "🧬 Variants after filtering: $variant_count"

        bcftools query -f '%CHROM\t%POS\n' "$vcf_pop_out" > "$vcf_var_out"

        LC_ALL=C awk 'NR==FNR{a[$1 FS $2]; next} ($1 FS $2) in a' "$vcf_var_out" "$dbSNPv157_pos" > "$flt_var_out"

        bcftools view -R "$flt_var_out" -Oz -o "$vcf_flt_out" "$vcf_pop_out"

        bcftools index "$vcf_flt_out"

        echo "✅ Done: chr${chr} ${pop}"
    } &> "$log_file"
}
export -f process_chr_pop
export cwd vcf_dir sample_dir dbSNPv157_pos vcf_pop_dir

# ---- PARALLEL EXECUTION ---- #
echo "🚀 Starting parallel processing..."

chroms=$(echo {1..22} X)

parallel -j 8 process_chr_pop ::: $chroms ::: "${populations[@]}"




# ----MERGE VCF FILE CHR/POP ---- #
# Base directory and config
export base="/user/home/xd14188/repo/1kg_hg38"
export vcf_pop_dir="${base}/vcfPops"
export outdir="${base}/mergedPopVcf"
export populations=("AFR" "EUR" "EAS" "SAS" "AMR")
mkdir -p "$outdir"

merge_one_pop() {
  pop="$1"
  echo "🔄 Merging $pop..."
  vcf_files=()
  for chr in {1..22} X; do
      vcf="${vcf_pop_dir}/chr${chr}_${pop}_only_filtered_SNP.vcf.gz"
      [[ -f "$vcf" ]] && vcf_files+=("$vcf")
  done
  [[ ${#vcf_files[@]} -eq 0 ]] && echo "⚠️ No files for $pop" && return
  out="${outdir}/${pop}_merged.vcf.gz"
  bcftools concat -Oz -o "$out" "${vcf_files[@]}"
  bcftools index --csi "$out"
  echo "✅ Done: $out"
}
export -f merge_one_pop
parallel -j 5 merge_one_pop ::: "${populations[@]}" 2>&1 | tee "${outdir}/merge_all.log"


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
