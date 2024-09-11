### Histograms of beta values for NonCpg methylatio (together and sep as CHG and CHH)

setwd("/mnt/data1/Eilis/Projects/Rat_Neuron/methylation_extraction/AllMethylation")
load("CpG_PresentAllSamples.Rdata")


setwd("/mnt/data1/Eilis/Projects/Rat_Neuron/")

### filter to min read depth of 10
all<-all[which(all[,34] >= 10),]

### histogram of beta values per sample

pdf("AllMethylation/plots/SummaryStats/Histogram_BetaValues_PerSample_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
par(mfrow = c(2,4))
for(i in (seq(from = 1, to = 29, by = 4))){
hist(all[,i], main = paste(gsub("_Methylation", "", colnames(all)[i]), "CpG Methylation", sep = " : "), xlab = "Beta Value", ylab = "", axes = FALSE)
axis(1, cex = 0.8)
axis(2, las = 2, format(c(0,100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000), scientific = FALSE))
}
dev.off()


pdf("AllMethylation/plots/SummaryStats/Histogram_BetaValues_PerSample_CpGMethylation_minReadDepth10_LogScale.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
par(mfrow = c(2,4))
for(i in (seq(from = 1, to = 29, by = 4))){
myhist<-hist(all[,i], plot = FALSE)
barplot(myhist$counts, main = paste(gsub("_Methylation", "", colnames(all)[i]), "CpG Methylation", sep = " : "), xlab = "Beta Value", log = "y", axes = FALSE, space = c(0,0))
axis(1, myhist$breaks[seq(from = 1, to = 21, by = 2)], at = c(0:(length(myhist$breaks)-1))[seq(from = 1, to = 21, by = 2)], tick = FALSE, cex = 0.8)
axis(2, las = 2, labels = format(c(1,10000, 50000, 100000, 500000), scientific = FALSE), at = c(1,10000, 50000, 100000, 500000))
}
dev.off()

### calc average methylation per site
meth.mean<-apply(all[,seq(from = 1, to = 29, by = 4)], 1, mean)
meth.median<-apply(all[,seq(from = 1, to = 29, by = 4)], 1, median)
names(meth.mean)<-rownames(all)
names(meth.median)<-rownames(all)


locs<-unlist(strsplit(rownames(all), "_"))
locs<-matrix(locs, ncol = 2, byrow = TRUE)
rownames(locs)<-rownames(all)


pdf("AllMethylation/plots/SummaryStats/Boxplot_MeanBetaValues_PerChr_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
boxplot(list(meth.mean[which(locs[,1] == "chr1")], meth.mean[which(locs[,1] == "chr2")], meth.mean[which(locs[,1] == "chr3")], meth.mean[which(locs[,1] == "chr4")],
meth.mean[which(locs[,1] == "chr5")],meth.mean[which(locs[,1] == "chr6")],meth.mean[which(locs[,1] == "chr7")],meth.mean[which(locs[,1] == "chr8")],
meth.mean[which(locs[,1] == "chr9")],meth.mean[which(locs[,1] == "chr10")],meth.mean[which(locs[,1] == "chr11")],meth.mean[which(locs[,1] == "chr12")], 
meth.mean[which(locs[,1] == "chr13")],meth.mean[which(locs[,1] == "chr14")],meth.mean[which(locs[,1] == "chr15")],meth.mean[which(locs[,1] == "chr16")],
meth.mean[which(locs[,1] == "chr17")],meth.mean[which(locs[,1] == "chr18")],meth.mean[which(locs[,1] == "chr19")],meth.mean[which(locs[,1] == "chr20")],
meth.mean[which(locs[,1] == "chrX")],meth.mean[which(locs[,1] == "chrM")]), names = c(1:20, "X", "M"), ylab = "Mean Beta Value", main = "CpG Methylation")

dev.off()


pdf("AllMethylation/plots/SummaryStats/Boxplot_MedianBetaValues_PerChr_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
boxplot(list(meth.median[which(locs[,1] == "chr1")], meth.median[which(locs[,1] == "chr2")], meth.median[which(locs[,1] == "chr3")], meth.median[which(locs[,1] == "chr4")],
meth.median[which(locs[,1] == "chr5")],meth.median[which(locs[,1] == "chr6")],meth.median[which(locs[,1] == "chr7")],meth.median[which(locs[,1] == "chr8")],
meth.median[which(locs[,1] == "chr9")],meth.median[which(locs[,1] == "chr10")],meth.median[which(locs[,1] == "chr11")],meth.median[which(locs[,1] == "chr12")], 
meth.median[which(locs[,1] == "chr13")],meth.median[which(locs[,1] == "chr14")],meth.median[which(locs[,1] == "chr15")],meth.median[which(locs[,1] == "chr16")],
meth.median[which(locs[,1] == "chr17")],meth.median[which(locs[,1] == "chr18")],meth.median[which(locs[,1] == "chr19")],meth.median[which(locs[,1] == "chr20")],
meth.median[which(locs[,1] == "chrX")],meth.median[which(locs[,1] == "chrM")]), names = c(1:20, "X", "M"), ylab = "Median Beta Value", main = "NonCpG Methylation")

dev.off()
