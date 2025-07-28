# Emma Walker
# E.M.Walker@exeter.ac.uk
# 06/02/19

### Manhattan plot for Isabel

setwd("/mnt/data1/EmmaW/Isabel/manhattan/")

library(ggplot2)
library(stringr)
library(qqman)

#get data
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeStats.Rdata")
genotype_stats <- as.data.frame(RRBSGenotypeStat)
colnames(genotype_stats) <- c("F_value", "Degree_of_freedom", "Effect_size", "p_value", "FDR", "n_Total")

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

### qq plot
qq(genotype_stats$p_value)

### manhattan plot
# rotate x axis labels so they all fit on the plot
#par(las = 2)
#mkae plot
manhattan(genotype_stats, col = rainbow(22), p = "p_value",
        chrlabs = c(1:19, "X", "Y", "MT"), annotatePval = 0.01, annotateTop = T,
        ylim = c(0,30), cex.axis = 0.6)


### some useful notes (see help(manhattan) for more detail)

# to change position of lines:
# genomewideline = xxx
# suggestiveline = xxx

# to highlight certain SNPs
# create a vector of SNPS (e.g. mySNPs), then:
#> highlight = mySNPs
