dat <- read.table("./sampleInfo.txt", header = T, sep = "\t")

print(unique(dat$Superpopulation.code))

dat[which(dat$Superpopulation.code == "EUR,AFR"), ]

data <- dat[-which(dat$Superpopulation.code == "EUR,AFR"), ]

print(unique(data$Superpopulation.code))
pops <- unique(data$Superpopulation.code)


for (i in 1:length(pops)) {
  out_name <- data[which(data$Superpopulation.code == pops[i]), "Sample.name"]
  print(length(out_name))
  # write.table(out_name, file = paste0("sampleName_", pops[i], ".txt"), quote = F, sep = "\n", row.names = F, col.names = F)
}