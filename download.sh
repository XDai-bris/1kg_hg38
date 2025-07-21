#!/bin/bash

#SBATCH --job-name=download_job
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:0:0
#SBATCH --mem=4G
#SBATCH --account=smed001801


download_http="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr"
gz_file_end=".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
tbi_file_end=".filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
vcf_file_dir="/user/work/xd14188/repo/1kg_hg38/vcfFiles"

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/20220804_manifest.txt

mkdir -p "${vcf_file_dir}"

for chr in 3 4 X; do
    wget "${download_http}${chr}${gz_file_end}" -P "${vcf_file_dir}"
    wget "${download_http}${chr}${tbi_file_end}" -P "${vcf_file_dir}"
done

# # Verify checksums
# grep -E "chr(3|4|X).*vcf.gz" 20220804_manifest.txt > manifest_subset.txt

# cd "${vcf_file_dir}"
# while read -r checksum filename; do
#     echo "$checksum  $filename" | md5sum -c -
# done < ../manifest_subset.txt
