# 1000 Genomes GRCh38 Sample Processing

## Overview

This project processes the `*.fam` from the [1000 Genomes GRCh38 data portal] to generate population-specific sample lists.

## Steps

### 1. Download Sample Metadata

Check if all 2,054 unrelated samples are in the latest 3,202 vcf files.

Download the `2054.txt` file containing metadata for **2,054 unrelated samples**:
📁 [Download here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index)

Download vcf files for **3,202 samples**
📁 [Download here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/*)
---

### 2. Check with old Hg37 version Fam
Download file (http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)
And check if the previous fam sample ID in the latest 3,202 vcf sample list.
```r
sam_3202 <- read.table("./3202_vcf_samples.txt", header = F)[, 1]

# List of populations
sam_2054_pop <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Loop through each population
for (pop in sam_2054_pop) {
  # Construct the input file path
  file_path <- paste0("./oldHg37VersionFim/", pop, ".fam")
  
  # Read the .fam file
  tmp <- read.table(file_path, header = FALSE)
  
  # Extract the sample IDs (second column)
  tmp_smp <- tmp$V2
  
  # Check if each sample is in sam_3202
  check <- tmp_smp %in% sam_3202
  
  # Print population and whether all samples matched
  cat("\nPopulation:", pop, " - all matched:", all(check), "\n")
  
  # If all matched, write tmp_smp to ./<pop>.fam
  if (all(check)) {
    out_path <- paste0("./", pop, ".fam")
    write.table(tmp_smp, file = out_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("Written to", out_path, "\n")
  }
}
```
# Console output block (triple backticks + text)
Population: AFR  - all matched: TRUE 
Written to ./AFR.fam 

Population: AMR  - all matched: TRUE 
Written to ./AMR.fam 

Population: EAS  - all matched: TRUE 
Written to ./EAS.fam 

Population: EUR  - all matched: TRUE 
Written to ./EUR.fam 

Population: SAS  - all matched: TRUE 
Written to ./SAS.fam 

- `AFR` — African Ancestry  
- `AMR` — Admixed American  
- `EAS` — East Asian  
- `EUR` — European  
- `SAS` — South Asian

