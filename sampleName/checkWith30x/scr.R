vcf_x30_name <- read.table("./sampleInfo_30x_vcf.txt", header = F)[, 1]

info_x20 <- read.table("../sampleInfo.txt", header = T, sep = "\t")
info_x20_name <- info_x20$Sample.name
info_x20_bioID <- info_x20$Biosample.ID

info_x30 <- read.table("./sampleInfo_30x.txt", header = T, sep = "\t")
info_x30_name <- info_x30$Sample.name
info_x30_bioID <-  info_x30$Biosample.ID



out_x20_From_x30 <- info_x20_name[!info_x20_name %in% info_x30_name]; print(out_x20_From_x30) 
write.table(out_x20_From_x30, file = "./x20_sample_NOT_in_x30.txt", row.names = F, col.names = F, quote = F)


out_x20_From_x30_bioID <- info_x20_bioID[!info_x20_bioID %in% info_x30_bioID]; print(out_x20_From_x30_bioID)
write.table(out_x20_From_x30_bioID, file = "./x20_BioID_NOT_in_x30.txt", row.names = F, col.names = F, quote = F)

out_x30_From_x30vcf <- vcf_x30_name[!info_x30_name %in% vcf_x30_name]; out_x30_From_x30vcf
# character(0)
all(out_x30_From_x30vcf); length(vcf_x30_name) == length(vcf_x30_name)
# [1] TRUE
# [1] TRUE