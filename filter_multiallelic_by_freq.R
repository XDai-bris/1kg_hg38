#!/usr/bin/env Rscript

# === Parse Arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript filter_multiallelic_by_altfreq.R <bim_file> <afreq_file> <output_prefix>")
}

bim_file <- args[1]
afreq_file <- args[2]
output_prefix <- args[3]

# === Libraries ===
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# === Step 1: Read AFREQ ===
message("[1/5] Reading allele frequency file: ", afreq_file)
afreq <- read_tsv(afreq_file, col_names = FALSE, comment = "#", col_types = cols(
  X1 = col_character(),   # CHR
  X2 = col_character(),   # ID
  X3 = col_character(),   # REF
  X4 = col_character(),   # ALT
  X5 = col_character(),   # PROVISIONAL_REF
  X6 = col_double(),      # ALT_FREQ
  X7 = col_integer()      # OBS_CT
))
colnames(afreq) <- c("CHR", "ID", "REF", "ALT", "PROVISIONAL_REF", "ALT_FREQ", "OBS_CT")

# === Step 2: Read BIM ===
message("[2/5] Reading BIM file: ", bim_file)
bim <- read_table(bim_file,
  col_names = c("CHR", "ID", "CM", "POS", "ALT", "REF"),
  col_types = cols(
    CHR = col_character(),
    ID  = col_character(),
    CM  = col_double(),
    POS = col_integer(),
    ALT = col_character(),
    REF = col_character()
  )
)

# === Step 3: Merge by ID ===
message("[3/5] Merging BIM with AFREQ by variant ID...")
merged <- inner_join(bim, afreq, by = "ID") %>%
  rename(
    CHR = CHR.x,
    ALT = ALT.x,
    REF = REF.x
  )

# === Step 4: Keep highest ALT_FREQ per CHR:POS ===
message("[4/5] Selecting max ALT_FREQ per CHR:POS...")
filtered <- merged %>%
  filter(!is.na(ALT_FREQ)) %>%
  group_by(CHR, POS) %>%
  slice_max(order_by = ALT_FREQ, n = 1, with_ties = FALSE) %>%
  ungroup()

# === Step 5: Output files ===
keep_ids <- filtered$ID
remove_ids <- setdiff(merged$ID, keep_ids)

message("[5/5] Writing outputs...")
write_lines(keep_ids, paste0(output_prefix, "_keep.txt"))
write_lines(remove_ids, paste0(output_prefix, "_remove.txt"))

write_tsv(
  filtered %>% select(CHR, ID, CM, POS, ALT, REF),
  paste0(output_prefix, "_filtered.bim"),
  col_names = FALSE
)

message("✅ Done. Variants kept: ", length(keep_ids), " | removed: ", length(remove_ids))
