# 1000 Genomes GRCh38 Sample Processing

## Overview

This project processes the `sampleInfo.txt` from the [1000 Genomes GRCh38 data portal](https://www.internationalgenome.org/data-portal/data-collection/grch38) to generate population-specific sample lists.

## Steps

### 1. Download Sample Metadata

Download the `sampleInfo.txt` file containing metadata for **2,709 samples**:

📁 [Download here](https://www.internationalgenome.org/data-portal/data-collection/grch38)

---

### 2. Group Samples by Superpopulation

Using the `Superpopulation.code` column, group samples into the following five superpopulations:

- `AFR` — African Ancestry  
- `AMR` — Admixed American  
- `EAS` — East Asian  
- `EUR` — European  
- `SAS` — South Asian

---

### 3. Filter Multi-Population Sample

One sample was found with multiple superpopulation labels:
```r
> unique(dat$Superpopulation.code)
[1] "EUR" "EAS" "AMR" "AFR" "SAS" "EUR,AFR"
> dat[which(dat$Superpopulation.code == "EUR,AFR"), ]
  Sample.name   Sex  Biosample.ID  Population.code  Superpopulation.code
  HG01783       male SAME124427    IBS,MSL          EUR,AFR
```


### 4. Output Files

After filtering out the multi-superpopulation sample, the following files are generated:

#### Population-specific sample lists:

- `sampleName_AFR.txt`
- `sampleName_AMR.txt`
- `sampleName_EAS.txt`
- `sampleName_EUR.txt`
- `sampleName_SAS.txt`

#### Combined sample list:

- `sampleName.txt` — contains all filtered samples

Each file contains one `Sample.name` per line for the corresponding superpopulation.


# Below: Process to Check Sample Omitted from 20x Dataset to 30x Dataset
# Sample Information Processing for 1000 Genomes 30x GRCh38 Dataset (in R)

This repository documents the steps to verify sample inclusion across datasets based on the **30x high-coverage GRCh38 phased data** from the [International Genome Sample Resource (IGSR)](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

All processing is done using **R**.

---

## Steps

### 1. Download and Load the 30x Sample Info List

Download the sample info metadata file from the 30x GRCh38 data portal and save it as `sampleInfo_30x.txt`.

Load 30x sample list:

```r
info_x30 <- read.table("./sampleInfo_30x.txt", header = T, sep = "\t")
info_x30_name <- info_x30$Sample.name
info_x30_bioID <-  info_x30$Biosample.ID
```

---

### 2. Check if Local Samples Exist in the 30x Sample List

Load your local sample list:

```r
info_x20 <- read.table("../sampleInfo.txt", header = T, sep = "\t")
info_x20_name <- info_x20$Sample.name
info_x20_bioID <- info_x20$Biosample.ID
```

---

### 3. Extract Sample Names from VCF File Using `bcftools`

Run this **outside R** in the shell to extract sample names from the VCF:

```bash
bcftools query -l ../../vcfFiles/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz > sampleInfo_30x_vcf.txt
```

Then load the list in R:

```r
vcf_x30_name <- read.table("./sampleInfo_30x_vcf.txt", header = F)[, 1]
```

---

### 4. Check 30x Sample List vs 30x VCF Sample Names

Check if all samples in `sampleInfo_30x` are in the VCF file `sampleInfo_30x_vcf.txt`:

```r
out_x30_From_x30vcf <- vcf_x30_name[!info_x30_name %in% vcf_x30_name]; out_x30_From_x30vcf
# character(0)
all(out_x30_From_x30vcf); length(vcf_x30_name) == length(vcf_x30_name)
# [1] TRUE
# [1] TRUE
```

---

### 5. Check Local Samples vs VCF Sample Names

Check if local samples are in the VCF list:

```r
out_x20_From_x30 <- info_x20_name[!info_x20_name %in% info_x30_name]; print(out_x20_From_x30) 
write.table(out_x20_From_x30, file = "./x20_sample_NOT_in_x30.txt", row.names = F, col.names = F, quote = F)


out_x20_From_x30_bioID <- info_x20_bioID[!info_x20_bioID %in% info_x30_bioID]; print(out_x20_From_x30_bioID)
write.table(out_x20_From_x30_bioID, file = "./x20_BioID_NOT_in_x30.txt", row.names = F, col.names = F, quote = F)
```
Threre are 137 samples in 20x are NOT in 30x!
Saved those 137 samples' sample ID and Bio ID in files `x20_sample_NOT_in_x30.txt` and `x20_sample_NOT_in_x30.txt`
---
