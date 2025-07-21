dat <- read.table("./sampleInfo.txt", header = T, sep = "\t")

print(unique(dat$Superpopulation.code))

dat[which(dat$Superpopulation.code == "EUR,AFR"), ]

data <- dat[-which(dat$Superpopulation.code == "EUR,AFR"), ]

# > print(unique(dat$Superpopulation.code))
# [1] "EUR"     "EAS"     "AMR"     "AFR"     "SAS"     "EUR,AFR"
# > 
#   > dat[which(dat$Superpopulation.code == "EUR,AFR"), ]
# Sample.name  Sex Biosample.ID Population.code Population.name Superpopulation.code               Superpopulation.name
# 1785     HG01783 male   SAME124427         IBS,MSL   Iberian,Mende              EUR,AFR European Ancestry,African Ancestry
# Population.elastic.ID                                                               Data.collections
# 1785               IBS,MSL 1000 Genomes on GRCh38,1000 Genomes 30x on GRCh38,1000 Genomes phase 3 release

print(unique(data$Superpopulation.code))
pops <- unique(data$Superpopulation.code)


for (i in 1:length(pops)) {
  out_name <- data[which(data$Superpopulation.code == pops[i]), "Sample.name"]
  print(length(out_name))
  # write.table(out_name, file = paste0("sampleName_", pops[i], ".txt"), quote = F, sep = "\n", row.names = F, col.names = F)
}