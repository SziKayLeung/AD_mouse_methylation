# Isabel Castanho / Emma Walker
# I.Castanho@exeter.ac.uk/ E.M.Walker@exeter.ac.uk
# 06/02/19

### Manhattan plots

setwd("/mnt/data1/isabel/RRBS/")

#library(ggplot2)
library(stringr)
library(qqman)
library(GenABEL)

color_Tg4510_TG <- "#00AEC9"
color_J20_TG <- "#FF5A62"

##### rTg4510 #####
## Genotype ##
# get data
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeStats.Rdata")
genotype_stats <- as.data.frame(RRBSGenotypeStat)

# add in columns needed to make manhattan plot
genotype_stats$SNP <- row.names(genotype_stats)
genotype_stats$CHR <- sapply(1:nrow(genotype_stats), function(i) str_extract(genotype_stats$SNP[i], "chr[^_]*"))
genotype_stats$CHR <- substring(genotype_stats$CHR, 4)
genotype_stats$BP <-sub(".*_", "", genotype_stats$SNP)

# change chrX to 20, chrY to 21 and M to 22
genotype_stats$CHR <- gsub("X", 20, genotype_stats$CHR)
genotype_stats$CHR <- gsub("Y", 21, genotype_stats$CHR)   
genotype_stats$CHR <- gsub("M", 22, genotype_stats$CHR)
genotype_stats$CHR <- as.numeric(genotype_stats$CHR)
genotype_stats$BP <- as.numeric(genotype_stats$BP)

# remove genes disrupted due to rTg4510 transgenes
row_names_remove <- c("chr12_115283603", "chr12_115561115", "chr12_115626480", "chr12_115673830", "chr12_116078466", "chr12_116078645", "chr12_116281566",
	"chr12_116399076", "chr12_11640163", "chr12_116403137", "chr12_116404504", "chr12_116404515", "chr12_116404525", "chr12_116404531", "chr12_116404566",
	"chr12_116404588", "chr12_116405400", "chr12_116488901", "chr12_116641364", "chr12_116641377", "chr12_116767580", "chr12_116979647", "chr12_117155000",
	"chr12_117155978", "chr12_117169957", "chr12_117177808", "chr12_117177849", "chr12_117179931", "chr12_117196967", "chr12_117200460", "chr12_117230891",
	"chr12_117258087", "chr12_117268435", "chr2_131909895", "chr2_131909902", "chr2_131909942", "chr2_131909943", "chr2_131909951", "chr2_131909953",
	"chr2_131909959", "chr2_131909972", "chr2_131909982", "chr2_131909987", "chr2_131909990", "chr2_131910163", "chr2_131910165", "chr2_131910181")

genotype_stats_clean <- genotype_stats[!(row.names(genotype_stats) %in% row_names_remove), ]


# Saving on object in RData format
#save(genotype_stats, genotype_stats_clean, file = "/mnt/data1/isabel/RRBS/rTg4510GenotypeStatsManhattan.RData")
# To load the data again
#load("/mnt/data1/isabel/RRBS/rTg4510GenotypeStatsManhattan.RData")

nrow(genotype_stats) # 1066467 
nrow(genotype_stats_clean) # 1066420

### qq plot
png(file = "/mnt/data1/isabel/RRBS/rTg4510_genotype_qqplot_clean.png", width=500, height=500)
qq(genotype_stats_clean$FDR_p)
dev.off()
# Genomic inflation factor lambda
data <- genotype_stats_clean$FDR_p
estlambda(data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1) # estimate = 1.035112

png(file = "/mnt/data1/isabel/RRBS/rTg4510_genotype_qqplot_full.png", width=500, height=500)
qq(genotype_stats$FDR_p)
dev.off()
# calculate lambda
data <- genotype_stats$FDR_p
estlambda(data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1) # estimate = 1.045509


### manhattan plot filtered
# make plot
threshold <- 0.05
png(file = "/mnt/data1/isabel/RRBS/rTg4510_genotype_manhattan_clean.png",
	width=1500, height=500)
manhattan(genotype_stats_clean,
	p = "FDR_p",
	ylim = c(0,20.5),
	cex = 0.8, # point size
	cex.axis = 1.5,
	cex.lab = 1.2,
	col = c("black", color_Tg4510_TG), #col = rainbow(22),
	suggestiveline = FALSE,
	genomewideline = -log10(threshold),
	chrlabs = c(1:19, "X", "Y", "MT"),
	#annotatePval = 0.001,
	las = 2, # rotate x axis labels so they all fit on the plot
	annotateTop = FALSE)
dev.off()

### manhattan plot unfiltered
# make plot
threshold <- 0.05
png(file = "/mnt/data1/isabel/RRBS/rTg4510_genotype_manhattan_full.png",
	width=1500, height=500)
manhattan(genotype_stats,
	p = "FDR_p",
	ylim = c(0,30),
	cex = 0.8, # point size
	cex.axis = 1.5,
	cex.lab = 1.2,
	col = c("black", color_Tg4510_TG), #col = rainbow(22),
	suggestiveline = FALSE,
	genomewideline = -log10(threshold),
	chrlabs = c(1:19, "X", "Y", "MT"),
	#annotatePval = 0.001,
	las = 2, # rotate x axis labels so they all fit on the plot
	annotateTop = FALSE)
dev.off()



##### J20 #####
## Genotype ##
#get data
load(file="/mnt/data1/Thea/IsabelSamples/data/J20GenotypeStats.Rdata")
genotype_stats <- as.data.frame(RRBSGenotypeStat)

# add in columns needed to make manhattan plot
genotype_stats$SNP <- row.names(genotype_stats)
genotype_stats$CHR <- sapply(1:nrow(genotype_stats), function(i) str_extract(genotype_stats$SNP[i], "chr[^_]*"))
genotype_stats$CHR <- substring(genotype_stats$CHR, 4)
genotype_stats$BP <-sub(".*_", "", genotype_stats$SNP)

# change chrX to 20, chrY to 21 and M to 22
genotype_stats$CHR <- gsub("X", 20, genotype_stats$CHR)
genotype_stats$CHR <- gsub("Y", 21, genotype_stats$CHR)   
genotype_stats$CHR <- gsub("M", 22, genotype_stats$CHR)
genotype_stats$CHR <- as.numeric(genotype_stats$CHR)
genotype_stats$BP <- as.numeric(genotype_stats$BP)

# Saving on object in RData format
save(genotype_stats, file = "/mnt/data1/isabel/RRBS/J20GenotypeStatsManhattan.RData")
# To load the data again
#load("/mnt/data1/isabel/RRBS/J20GenotypeStatsManhattan.RData")

nrow(genotype_stats) # 

### qq plot
png(file = "/mnt/data1/isabel/RRBS/J20_genotype_qqplot.png", width=500, height=500)
qq(genotype_stats$FDR_p)
dev.off()
# calculate lambda
data <- genotype_stats$FDR_p
estlambda(data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1) # estimate = 0.9198793


### manhattan plot
# make plot
threshold <- 0.05
png(file = "/mnt/data1/isabel/RRBS/J20_genotype_manhattan.png",
	width=1500, height=500)
manhattan(genotype_stats,
	p = "FDR_p",
	ylim = c(0,10),
	cex = 0.8, # point size
	cex.axis = 1.5,
	cex.lab = 1.2,
	col = c("black", color_J20_TG), # col = rainbow(22),
	suggestiveline = FALSE,
	genomewideline = -log10(threshold),
	chrlabs = c(1:19, "X", "Y", "MT"),
	#annotatePval = 0.001,
	las = 2, # rotate x axis labels so they all fit on the plot
	annotateTop = FALSE)
dev.off()


### some useful notes (see help(manhattan) for more detail)
# to change position of lines:
# genomewideline = xxx
# suggestiveline = xxx

# to highlight certain SNPs
# create a vector of SNPS (e.g. mySNPs), then:
#> highlight = mySNPs
