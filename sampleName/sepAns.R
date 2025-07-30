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
