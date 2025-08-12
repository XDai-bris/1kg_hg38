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
# Read with the file's own headers; keep '#CHROM' then normalize
afreq <- readr::read_tsv(afreq_file, show_col_types = FALSE)

# Normalize headers from PLINK2:
if ("#CHROM" %in% names(afreq)) names(afreq)[names(afreq) == "#CHROM"] <- "CHR"
if ("ALT_FREQS" %in% names(afreq)) names(afreq)[names(afreq) == "ALT_FREQS"] <- "ALT_FREQ"
if (!"PROVISIONAL_REF" %in% names(afreq)) afreq$PROVISIONAL_REF <- NA_character_

# Enforce types
afreq <- dplyr::mutate(
  afreq,
  CHR = as.character(CHR),
  ID = as.character(ID),
  REF = as.character(REF),
  ALT = as.character(ALT),
  PROVISIONAL_REF = as.character(PROVISIONAL_REF),
  ALT_FREQ = as.numeric(ALT_FREQ),
  OBS_CT = as.integer(OBS_CT)
)

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
afreq_dedup <- afreq %>%
  group_by(ID) %>%
  slice_max(order_by = ALT_FREQ, n = 1, with_ties = FALSE) %>%
  ungroup()

merged <- inner_join(bim, afreq_dedup, by = "ID") %>%
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
