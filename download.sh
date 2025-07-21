

#!/bin/bash

conda activate baseToolsStation

download_http="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr"
gz_file_end=".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
tbi_file_end=".filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
vcf_file_dir="/Users/xd14188/Desktop/UoB/1kg_hg38/vcfFiles"

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/20220804_manifest.txt

mkdir -p "${vcf_file_dir}"

for chr in 1..22; do
    wget "${download_http}${chr}${gz_file_end}" -P "${vcf_file_dir}"
    wget "${download_http}${chr}${tbi_file_end}" -P "${vcf_file_dir}"
done

wget "${download_http}X.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" -P "${vcf_file_dir}"
wget "${download_http}X.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi" -P "${vcf_file_dir}"

./vcfFiles/check_md5.sh
