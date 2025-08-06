# 1KGP High-Coverage VCF Processing Pipeline (hg38, HPC-Optimized)

This HPC-ready pipeline automates extraction, filtering, and transformation of high-coverage 1000 Genomes Project (hg38) VCFs into **population-specific**, **rsID-verified**, and **PLINK-compatible** files.\
Built for high-performance computing environments using `plink2`, `bcftools`, R, and `GNU parallel`.

---

## 📁 Project Structure

```
1kg_hg38/
├── vcfFiles/                  # Input: VCF files per chromosome (chr1–22, X)
├── sampleName/               # Input: Sample list per population (e.g., AFR.fam)
├── bedFiles_maf001/          # Output: Final cleaned PLINK files (per chr × pop)
│   └── chr<chr>_<pop>/       # → final_output.{bed,bim,fam}
├── tmp/                      # Intermediate files per job
├── logs/                     # Per-task logs from chr/pop processing
├── filterVar/                # Variants that couldn't be renamed to rsIDs
├── fltMultiallelics/         # Multiallelic variants filtered out
├── fltUnMatchedRsID/         # Clean variants with no rsID match
├── scripts/                  # Optional: scripts used in the process
├── mergeChrBed.sh            # Script to merge across chromosomes per pop
├── filter_multiallelic_by_freq.R  # R script to filter Multiallelics by ALT_FREQ
└── run_chr_pop_pipeline_imputedVarID..sh   # Main pipeline driver for Chr1-22 ChrX in HPC
```

---

## 🧰 Dependencies

Make sure the following tools are installed and available:

- [`plink2`](https://www.cog-genomics.org/plink/2.0/)
- [`bcftools`](https://samtools.github.io/bcftools/)
- [`GNU parallel`](https://www.gnu.org/software/parallel/)
- `R` with packages: `dplyr`, `readr`

---

## 🦬 Inputs

- **Chromosome VCFs**\
  One per chromosome (e.g., `chr10.vcf.gz`), downloaded from:

  ```
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
  ```

- **Sample Lists**\
  One `.fam` file per population, e.g.:

  ```
  sampleName/AFR.fam
  sampleName/EUR.fam
  ...
  ```

- **dbSNP rsID Map (per chromosome)**\
  Files like `chr21_dbsnp.tsv`, in format:

  ```
  CHR:POS:REF:ALT    rsID
  ```

  Located in:

  ```
  /user/home/xd14188/repo/data/genomeRef/cpra_rsID/
  ```

---

## ⚙️ Pipeline Workflow

### 1. **Per-Chromosome, Per-Population Filtering**

For each chromosome × population:

- Filters to population-specific individuals
- Impute mismatched/missing variant IDs in BIM file as CHR:POS:REF:ALT
- Filters variants to keep variants with **MAF ≥ 0.01**
- Computes alt allele frequency
- Filters out multiallelic sites using Max ALT frequency selection in R `filter_multiallelic_by_freq.R`
- Renames variants to `rsIDs` using exact CHR\:POS\:REF\:ALT match per CHR 
- Saves unmatched variants and removed multiallelics separately
- Final output: `final_output.{bed,bim,fam}` with rsIDs
- Verified `final_output.{bed,bim,fam}` via allele frequency (`--freq`) calculation per merged output

🌟 Output directory:\
`bedFiles_maf001/chr<chr>_<pop>/final_output.*`

---

### 2. **Parallelized Execution**

Run across multiple chromosome × population pairs via:

> Uses GNU `parallel` to dispatch jobs concurrently. Modify `-j 12` in the script to suit your HPC's limits.

---

### 3. **Merging Across Chromosomes (Optional)**

Merge PLINK files across all chromosomes per population using:

```bash
./mergeChrBed.sh
```

Output saved to `merged_genome/<POP>_merged.*`

Also verifies merge via allele frequency (`--freq`) per merged output.

---

## 🧪 Validation

After processing, each final PLINK file is validated using `plink2 --freq` to ensure integrity.

Each folder contains:

- `final_output.bed/bim/fam`
- `final_output_freq.afreq` (for quality check)
- Log: `final_output.log`

---

## 🌍 Supported Populations

- `EUR` – European
- `AFR` – African
- `EAS` – East Asian
- `SAS` – South Asian
- `AMR` – Admixed American

---

## 📄 Key Outputs

| Folder              | Contents                                         |
| ------------------- | ------------------------------------------------ |
| `bedFiles_maf001/`  | Final PLINK files per chr/pop (`final_output.*`) |
| `fltMultiallelics/` | BIM files with removed multiallelic variants     |
| `fltUnMatchedRsID/` | BIM entries that could not be matched to rsID    |
| `filterVar/`        | All unmatched entries not kept in final output   |
| `merged_genome/`    | Per-population merged genome-wide PLINK files    |

---


## 👎 Notes

- 1000G vcfs differs from the common multiallelic representation (e.g., using commas ',' to separate ALT alleles), which can be easily overlooked during processing:
  - SNVs and INDELs were phased in SHAPEIT2-duohmm, which does not handle multiallelic variant
  - we first split `multiallelic variants into separate rows` while left-aligning and normalizing INDELs using bcftools norm (Li, 2011). To phase both biallelic and multiallelic variants, we `shifted the positions of additional ALT alleles (2nd, 3rd, etc.) by 1 or more base pairs` to ensure unique start positions required by SHAPEIT2. After phasing, `positions were shifted back to their original coordinates`.
- Using `awk` to reading chr_bim and check with c:p:r:a format is per chr is much faster
- ChrX without variant ID in the VCF, which needs make c:p:r:a ID at first BIM file
- Logs are stored in `logs/` with full stdout/stderr per job
- You can re-run individual jobs easily by modifying `run_chr_pop_pipeline.sh`

---

## 📄 Citation

If this pipeline helps your work, please cite:

- The [1000 Genomes Project](https://www.internationalgenome.org/)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [bcftools](http://samtools.github.io/bcftools/)
- [GNU Parallel](https://www.gnu.org/software/parallel/)

---

## 👨‍💻 Maintainer

**Xiaoyang Dai**\
📧 [xd14188@bristol.ac.uk](mailto\:xd14188@bristol.ac.uk)

