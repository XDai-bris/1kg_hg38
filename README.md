# 1KGP High-Coverage VCF Sample Processor (hg38)

This script is designed to automate the extraction and filtering of 1000 Genomes Project high-coverage VCF files (hg38) by population, using bcftools and plink. It supports sample-specific filtering based on population and applies dbSNP filters to produce PLINK-compatible BED files for downstream genomic analysis.

---

## 📁 Project Structure

```
tools/genomeRef/            # Reference files (e.g., dbSNPv157_hg38)
1kg_hg38/
├── vcfFiles/               # Input: VCF files (chr1..22)
├── sampleName/             # Input: Sample lists per population (e.g., sampleName_EUR.txt)
├── tmp_vcf_sampleName/     # Temp: Sample names extracted from VCFs
├── notInVcfSample/         # Output: Samples missing from VCFs
├── bedFiles_maf001/        # Output: Population-filtered BED files with MAF > 0.01
├── vcfPops/                # Output: Population-filtered VCFs
```

---

## 🧰 Dependencies

Make sure the following tools are installed and available in your environment:

- [`bcftools`](https://samtools.github.io/bcftools/)
- [`plink`](https://www.cog-genomics.org/plink/), version PLINK1.9 as PLINK2 can not be installed
- [`GNU parallel`](https://www.gnu.org/software/parallel/)
- [`20220422_3202_phased_SNV_INDEL_SV`](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
---

## 🗂 Inputs

- **VCF files:**  
  From the 1000 Genomes Project high-coverage data [`20220422_3202_phased_SNV_INDEL_SV`] for each chromosomes 1–22,X located in `vcfFiles/`.

- **Sample Lists:**  
  One per population (e.g., `sampleName_EUR.txt`, `sampleName_AFR.txt`) stored in `sampleName/`.
    this should be prepared, details see README.md in `sampleName/`
- **dbSNP v157 hg38 BED file:**  
  Located at `tools/genomeRef/dbSNPv157_hg38/dbsnp157_biallelic_rs_hg38.bed`. Used for SNP filtering.
    this should be prepared, details see README.md in `tools/genomeRef/dbSNPv157_hg38/`
---

## ⚙️ What It Does

1. **Extracts sample names** from VCFs for each chromosome (1–22, X).
2. **Verifies consistency** of sample names across each chromosomes vcf.gz files(using chr1 as reference).
3. **Filters VCF files**:
   - Keeps only samples for a given population.
   - Filters variants based on dbSNP BED file.
   - Outputs compressed VCFs and PLINK BED files with MAF > 0.01.
4. **Parallel processing** of chromosomes and populations (via GNU `parallel`).

---

## 🚀 Usage
    select suitable parallel -j 3 at end of `plink_parallel.sh`
    chmod +x plink_parallel.sh
    ./plink_parallel.sh

---

## 🌍 Supported Populations

- EUR – European
- AFR – African
- EAS – East Asian
- SAS – South Asian
- AMR – Admixed American

---

## 📌 Notes

- The script currently only processes **chr1** for testing. To process all chromosomes, update:
  ```bash
  parallel -j 3 process_chr_pop ::: {1..22} ::: EUR AFR EAS SAS AMR
  ```
- Make sure your `sampleName_*.txt` files contain sample IDs matching those in the VCF files.
- There are some samples are not in Vcf but in the different sampleName_Pop
---

## 📤 Output

- `vcfPops/chr<chr>_<pop>_only.vcf.gz` – population-specific, dbSNP-filtered VCFs.
- `bedFiles_maf001/chr<chr>_<pop>_MAF01.*` – PLINK-formatted genotype data for MAF > 0.01.
- `notInVcfSample/*` - Samples in different sampleName_Pop but can not be found in chr1..22,X Vcf file in `vcfFiles/`
---


