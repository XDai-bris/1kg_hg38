# 1KGP High-Coverage VCF Processing Pipeline (hg38, HPC Version)

This HPC-optimized pipeline automates the extraction, filtering, and transformation of 1000 Genomes Project high-coverage (hg38) VCFs into population-specific, dbSNP-verified, and PLINK-compatible files. It leverages `bcftools`, `plink`, and GNU `parallel` for efficient chromosome- and population-level processing on high-performance clusters.

---

## 📁 Project Structure

```
1kg_hg38/
├── vcfFiles/              # Input: Downloaded Chromosome VCFs (chr1..22, X)
├── sampleName/            # Input: Sample list per population (e.g., sampleName_EUR.txt)
├── tmp_vcf_sampleName/    # Temp: Extracted sample names from input VCFs
├── notInVcfSample/        # Output: Samples listed but not found in VCFs
├── vcfPops/               # Output: Population-filtered, dbSNP-verified VCFs
│   └── logs/              # Logs for each chr/pop processing step
├── mergedPopVcf/          # Output: Merged population-level VCFs
└── mergedPopBed/          # Output: Final PLINK .bed/.bim/.fam files per population
```

---

## 🧰 Dependencies

Ensure the following tools are available in your environment:

- [`bcftools`](bcftools/1.19-openblas-2333)
- [`plink2`](plink2/2.00a4.3-netlib-lapack-qbk6) 
- [`GNU parallel`](http://www.gnu.org/software/parallel)
- High-coverage 1000 Genomes data:  
  [`20220422_3202_phased_SNV_INDEL_SV`](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)

---

## 🗂 Inputs

- **Chromosome VCFs:**  
  Used /script/download.sh to download from the 1000 Genomes high-coverage data (hg38), one file per chr1–22 and X, stored in `vcfFiles/`. 
- **Sample Lists:**  
  One per population (e.g., `sampleName_EUR.txt`, `sampleName_AFR.txt`), stored in `sampleName/`.  
  ➤ *Each list contains sample IDs corresponding to individuals from that population.*

- **dbSNP v157 hg38 (biallelic rsIDs only):**  
  TAB-delimited file: `dbsnp157_biallelic_rs_hg38_UCSC.tab`  
  Path: `repo/data/genomeRef/dbsnp157_biallelic_rs_hg38_UCSC.tab/`  
  ➤ *Used to retain known and reliable SNPs.*

---

## ⚙️ Pipeline Overview

### Step 1: **Chromosome + Population Filtering**
- Extracts VCF sample names
- Filters to keep only samples for the selected population
- Filters for SNPs only, and retains those with **MAF > 0.01**
- Outputs filtered `.vcf.gz` and associated variant list

### Step 2: **dbSNP Verification**
- Intersects filtered variants with `dbsnp157_biallelic_rs_hg38` to ensure known SNPs
- Generates:
  - Filtered `.vcf.gz`
  - `.tab` files listing positions before and after dbSNP filtering

### Step 3: **Parallel Execution**
- Steps 1–2 are run in parallel across all chromosomes and populations using `GNU parallel`

### Step 4: **VCF Merging**
- Merges filtered per-chromosome VCFs into population-wide VCFs
- Output: `mergedPopVcf/<pop>_merged.vcf.gz`

### Step 5: **Convert to PLINK Format**
- Uses `plink` to convert merged VCFs to `.bed/.bim/.fam` for downstream GWAS or analysis

---

## 🚀 Usage

1. **Make script executable:**
   ```bash
   chmod +x process_1kg_pipeline.sh
   ```

2. **Run the pipeline:**
   ```bash
   ./process_1kg_pipeline.sh
   ```

3. **Control parallelism:**  
   The number of concurrent jobs (`-j`) can be tuned at the bottom of the script to match your HPC cluster limits.

---

## 🌍 Supported Populations

- `EUR` – European  
- `AFR` – African  
- `EAS` – East Asian  
- `SAS` – South Asian  
- `AMR` – Admixed American  

---

## 📤 Output Files

- `vcfPops/chr<chr>_<pop>_only.vcf.gz` – SNPs filtered by population
- `vcfPops/chr<chr>_<pop>_only_varList.tab` – Variant positions (pre-dbsnp filtering)
- `vcfPops/chr<chr>_<pop>_only_varList_filtered.tab` – Variants retained after dbSNP filter
- `vcfPops/chr<chr>_<pop>_only_filtered_SNP.vcf.gz` – Final filtered VCF
- `mergedPopVcf/<pop>_merged.vcf.gz` – Merged across chromosomes
- `mergedPopBed/<pop>_plink.bed/.bim/.fam` – Final PLINK binary format
- `notInVcfSample/chr<chr>_notInVcfSample_<pop>.txt` – Missing samples

---

## ⚠️ Notes

- **Sample ID Matching:** Make sure your `sampleName_*.txt` contains correct sample IDs matching those in the 1KGP VCFs.
- **Missing Samples:** Some population lists may include samples not present in VCFs; those are logged in `notInVcfSample/`.
- **Logging:** Full per-task logs can be found in `vcfPops/logs/`.

---

## 📄 License & Citation

If using this pipeline for publication, cite the 1000 Genomes Project and tools (bcftools, plink and GUN parellel) as appropriate.

---
