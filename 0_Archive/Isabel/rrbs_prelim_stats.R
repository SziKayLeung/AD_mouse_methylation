#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS results table from Dorothea

# setwd
setwd("/mnt/data1/isabel/RRBS/")

############ rTg4510 ############
## Genotype rTg4510
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeStats.Rdata")
genotype_stats <- as.data.frame(RRBSGenotypeStat)
genotype_stats_FDR <- genotype_stats[which(genotype_stats[,"FDR_p"] < 0.05),]
write.csv(genotype_stats_FDR, file = "/mnt/data1/isabel/RRBS/rTg4510GenotypeStats_full.csv")

# remove genes disrupted due to rTg4510 transgenes
row_names_remove <- c("chr12_115283603", "chr12_115561115", "chr12_115626480", "chr12_115673830", "chr12_116078466", "chr12_116078645", "chr12_116281566",
	"chr12_116399076", "chr12_11640163", "chr12_116403137", "chr12_116404504", "chr12_116404515", "chr12_116404525", "chr12_116404531", "chr12_116404566",
	"chr12_116404588", "chr12_116405400", "chr12_116488901", "chr12_116641364", "chr12_116641377", "chr12_116767580", "chr12_116979647", "chr12_117155000",
	"chr12_117155978", "chr12_117169957", "chr12_117177808", "chr12_117177849", "chr12_117179931", "chr12_117196967", "chr12_117200460", "chr12_117230891",
	"chr12_117258087", "chr12_117268435", "chr2_131909895", "chr2_131909902", "chr2_131909942", "chr2_131909943", "chr2_131909951", "chr2_131909953",
	"chr2_131909959", "chr2_131909972", "chr2_131909982", "chr2_131909987", "chr2_131909990", "chr2_131910163", "chr2_131910165", "chr2_131910181")

genotype_stats_clean <- genotype_stats[!(row.names(genotype_stats) %in% row_names_remove), ]
genotype_stats_FDR_clean <- genotype_stats_clean[which(genotype_stats_clean[,"FDR_p"] < 0.05),]
write.csv(genotype_stats_FDR_clean, file = "/mnt/data1/isabel/RRBS/rTg4510GenotypeStats_clean.csv")


############ J20 ############
## Genotype J20
load(file="/mnt/data1/Thea/IsabelSamples/data/J20GenotypeStats.Rdata")
genotype_stats <- as.data.frame(RRBSGenotypeStat)
genotype_stats_FDR <- genotype_stats[which(genotype_stats[,"FDR_p"] < 0.05),]
write.csv(genotype_stats_FDR, file = "/mnt/data1/isabel/RRBS/J20GenotypeStats.csv")






#### IGNORE FOR NOW

## Genotype by age
files <- c("/mnt/data1/isabel/RRBS/rTg4510stats1.csv",
	"/mnt/data1/isabel/RRBS/rTg4510stats2.csv",
	"/mnt/data1/isabel/RRBS/rTg4510stats3.csv",
	"/mnt/data1/isabel/RRBS/rTg4510stats4.csv")

ages <- c("2", "4", "6", "8")

for (index in 1:length(ages)) {
	file <- files[index]
	age <- ages[index]
	data <- read.csv(file, stringsAsFactors=FALSE, row.names=1)
	colnames(data) <- c("F_value", "Degree_of_freedom", "Effect_size", "p_value", "FDR", "n_Total")
	stats_FDR <- data[which(data[,"FDR"] < 0.05),]
	write.csv(stats_FDR, file = paste("/mnt/data1/isabel/RRBS/rTg4510Stats_", age, "months.csv", sep=""))
}


## Genotype by age
files <- c("/mnt/data1/isabel/RRBS/J20stats1.csv",
	"/mnt/data1/isabel/RRBS/J20stats2.csv",
	"/mnt/data1/isabel/RRBS/J20stats3.csv",
	"/mnt/data1/isabel/RRBS/J20stats4.csv")

ages <- c("6", "8", "10", "12")

for (index in 1:length(ages)) {
	file <- files[index]
	age <- ages[index]
	data <- read.csv(file, stringsAsFactors=FALSE, row.names=1)
	colnames(data) <- c("t_value", "Degree_of_freedom", "Effect_size", "p_value", "FDR", "n_Total")
	stats_FDR <- data[which(data[,"FDR"] < 0.05),]
	write.csv(stats_FDR, file = paste("/mnt/data1/isabel/RRBS/J20Stats_", age, "months.csv", sep=""))
}



########### other info / data
# rTg4510GenotypeMatrixFiltered.Rdata
#load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeMatrixFiltered.Rdata")