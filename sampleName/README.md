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
