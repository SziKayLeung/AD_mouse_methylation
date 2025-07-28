### calculate sumary statistics separately for CpG and nonCpG methylation

setwd("/mnt/data1/Eilis/Projects/Rat_Neuron/methylation_extraction/AllMethylation")
### start with non CpG
cpg<-read.table("All_CpG_Sites_AcrossAllSamples.txt")
cpg<-paste(cpg[,2], cpg[,3], sep = "_")

##create list of probes for each file

CreateProbeList<-function(x){
  rownames(x)<-paste(x[,1], x[,2], sep = "_")
  return(x)
}


##Load Bedgraph files
A3<-read.table("A3.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A3<-CreateProbeList(A3)
A3<-A3[intersect(cpg, rownames(A3)),]
A4<-read.table("A4.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A4<-A4[intersect(cpg, rownames(A4)),]
A5<-read.table("A5.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A5<-CreateProbeList(A5)
A5<-A5[intersect(cpg, rownames(A5)),]
A6<-read.table("A6.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A6<-CreateProbeList(A6)
A6<-A6[intersect(cpg, rownames(A6)),]
A7<-read.table("A7.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A7<-CreateProbeList(A7)
A7<-A7[intersect(cpg, rownames(A7)),]
A8<-read.table("A8.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A8<-CreateProbeList(A8)
A8<-A8[intersect(cpg, rownames(A8)),]
A9<-read.table("A9.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A9<-CreateProbeList(A9)
A9<-A9[intersect(cpg, rownames(A9)),]
A10<-read.table("A10.new.1_val_1.fq.gz_bismark_pe.sam.bismark.cov", header = F)
A10<-CreateProbeList(A10)
A10<-A10[intersect(cpg, rownames(A10)),]

A3<-as.matrix(A3[,-c(1:3)])
A4<-as.matrix(A4[,-c(1:3)])
A5<-as.matrix(A5[,-c(1:3)])
A6<-as.matrix(A6[,-c(1:3)])
A7<-as.matrix(A7[,-c(1:3)])
A8<-as.matrix(A8[,-c(1:3)])
A9<-as.matrix(A9[,-c(1:3)])
A10<-as.matrix(A10[,-c(1:3)])

A3<-A3[intersect(cpg, rownames(A3)),]
A4<-A4[intersect(cpg, rownames(A4)),]
A5<-A5[intersect(cpg, rownames(A5)),]
A6<-A6[intersect(cpg, rownames(A6)),]
A7<-A7[intersect(cpg, rownames(A7)),]
A8<-A8[intersect(cpg, rownames(A8)),]
A9<-A9[intersect(cpg, rownames(A9)),]
A10<-A10[intersect(cpg, rownames(A10)),]

save(A3,A4,A5,A6,A7,A8,A9,A10, filename = "CpGMethylation.Rdata")

### calculate coverage:
CalcCoverage<-function(x, lab){
  x<-cbind(x, x[,2]+x[,3])
  colnames(x)<-paste(lab, c('Methylation', 'nMReads', 'nNMReads', 'TotReads'), sep = "_")
  return(x)
}


A3<-CalcCoverage(A3, lab = "A3")
A4<-CalcCoverage(A4, lab = "A4")
A5<-CalcCoverage(A5, lab = "A5")
A6<-CalcCoverage(A6, lab = "A6")
A7<-CalcCoverage(A7, lab = "A7")
A8<-CalcCoverage(A8, lab = "A8")
A9<-CalcCoverage(A9, lab = "A9")
A10<-CalcCoverage(A10, lab = "A10")

save(A3,A4,A5,A6,A7,A8,A9,A10, file = "CpGMethylation.Rdata")

### summary statistics per sample:

summary.id<-matrix(data = NA, nrow = 8, ncol = 11)
colnames(summary.id)<-c("Number of Sites", "Min Read Depth", "1st Quantile Read Depth", "Mean Read Depth", "Median Read Depth",  "3rd Quantile Read Depth", "Max Read Depth", "Percentage Sites:Read Depth 5",  "Percentage Sites:Read Depth 10", "Sites:Read Depth 5","Sites:Read Depth 10")
rownames(summary.id)<-c("A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")

summary.id[1,1]<-nrow(A3)
summary.id[2,1]<-nrow(A4)
summary.id[3,1]<-nrow(A5)
summary.id[4,1]<-nrow(A6)
summary.id[5,1]<-nrow(A7)
summary.id[6,1]<-nrow(A8)
summary.id[7,1]<-nrow(A9)
summary.id[8,1]<-nrow(A10)

summary.id[1,2]<-min(A3[,4])
summary.id[2,2]<-min(A4[,4])
summary.id[3,2]<-min(A5[,4])
summary.id[4,2]<-min(A6[,4])
summary.id[5,2]<-min(A7[,4])
summary.id[6,2]<-min(A8[,4])
summary.id[7,2]<-min(A9[,4])
summary.id[8,2]<-min(A10[,4])

summary.id[1,3]<-quantile(A3[,4], 0.25)
summary.id[2,3]<-quantile(A4[,4], 0.25)
summary.id[3,3]<-quantile(A5[,4], 0.25)
summary.id[4,3]<-quantile(A6[,4], 0.25)
summary.id[5,3]<-quantile(A7[,4], 0.25)
summary.id[6,3]<-quantile(A8[,4], 0.25)
summary.id[7,3]<-quantile(A9[,4], 0.25)
summary.id[8,3]<-quantile(A10[,4], 0.25)

summary.id[1,4]<-mean(A3[,4])
summary.id[2,4]<-mean(A4[,4])
summary.id[3,4]<-mean(A5[,4])
summary.id[4,4]<-mean(A6[,4])
summary.id[5,4]<-mean(A7[,4])
summary.id[6,4]<-mean(A8[,4])
summary.id[7,4]<-mean(A9[,4])
summary.id[8,4]<-mean(A10[,4])

summary.id[1,5]<-median(A3[,4])
summary.id[2,5]<-median(A4[,4])
summary.id[3,5]<-median(A5[,4])
summary.id[4,5]<-median(A6[,4])
summary.id[5,5]<-median(A7[,4])
summary.id[6,5]<-median(A8[,4])
summary.id[7,5]<-median(A9[,4])
summary.id[8,5]<-median(A10[,4])

summary.id[1,6]<-quantile(A3[,4], 0.75)
summary.id[2,6]<-quantile(A4[,4], 0.75)
summary.id[3,6]<-quantile(A5[,4], 0.75)
summary.id[4,6]<-quantile(A6[,4], 0.75)
summary.id[5,6]<-quantile(A7[,4], 0.75)
summary.id[6,6]<-quantile(A8[,4], 0.75)
summary.id[7,6]<-quantile(A9[,4], 0.75)
summary.id[8,6]<-quantile(A10[,4], 0.75)

summary.id[1,7]<-max(A3[,4])
summary.id[2,7]<-max(A4[,4])
summary.id[3,7]<-max(A5[,4])
summary.id[4,7]<-max(A6[,4])
summary.id[5,7]<-max(A7[,4])
summary.id[6,7]<-max(A8[,4])
summary.id[7,7]<-max(A9[,4])
summary.id[8,7]<-max(A10[,4])

summary.id[1,8]<-length(which(A3[,4] >= 5))/nrow(A3)*100
summary.id[2,8]<-length(which(A4[,4] >= 5))/nrow(A4)*100
summary.id[3,8]<-length(which(A5[,4] >= 5))/nrow(A5)*100
summary.id[4,8]<-length(which(A6[,4] >= 5))/nrow(A6)*100
summary.id[5,8]<-length(which(A7[,4] >= 5))/nrow(A7)*100
summary.id[6,8]<-length(which(A8[,4] >= 5))/nrow(A8)*100
summary.id[7,8]<-length(which(A9[,4] >= 5))/nrow(A9)*100
summary.id[8,8]<-length(which(A10[,4] >= 5))/nrow(A10)*100

summary.id[1,9]<-length(which(A3[,4] >= 10))/nrow(A3)*100
summary.id[2,9]<-length(which(A4[,4] >= 10))/nrow(A4)*100
summary.id[3,9]<-length(which(A5[,4] >= 10))/nrow(A5)*100
summary.id[4,9]<-length(which(A6[,4] >= 10))/nrow(A6)*100
summary.id[5,9]<-length(which(A7[,4] >= 10))/nrow(A7)*100
summary.id[6,9]<-length(which(A8[,4] >= 10))/nrow(A8)*100
summary.id[7,9]<-length(which(A9[,4] >= 10))/nrow(A9)*100
summary.id[8,9]<-length(which(A10[,4] >= 10))/nrow(A10)*100

summary.id[1,10]<-length(which(A3[,4] >= 5))
summary.id[2,10]<-length(which(A4[,4] >= 5))
summary.id[3,10]<-length(which(A5[,4] >= 5))
summary.id[4,10]<-length(which(A6[,4] >= 5))
summary.id[5,10]<-length(which(A7[,4] >= 5))
summary.id[6,10]<-length(which(A8[,4] >= 5))
summary.id[7,10]<-length(which(A9[,4] >= 5))
summary.id[8,10]<-length(which(A10[,4] >= 5))

summary.id[1,11]<-length(which(A3[,4] >= 10))
summary.id[2,11]<-length(which(A4[,4] >= 10))
summary.id[3,11]<-length(which(A5[,4] >= 10))
summary.id[4,11]<-length(which(A6[,4] >= 10))
summary.id[5,11]<-length(which(A7[,4] >= 10))
summary.id[6,11]<-length(which(A8[,4] >= 10))
summary.id[7,11]<-length(which(A9[,4] >= 10))
summary.id[8,11]<-length(which(A10[,4] >= 10))

write.csv(summary.id, "EachSampleReadDepthStatsAllSites_CpGsites.csv")

### Create List of probes covered in all samples
probes_all<-intersect(rownames(A3), rownames(A4))
probes_all<-intersect(probes_all, rownames(A5))
probes_all<-intersect(probes_all, rownames(A6))
probes_all<-intersect(probes_all, rownames(A7))
probes_all<-intersect(probes_all, rownames(A8))
probes_all<-intersect(probes_all, rownames(A9))
probes_all<-intersect(probes_all, rownames(A10))

##obtain total number of sites assessed across all samples
allposs<-unique(c(rownames(A3), rownames(A4), rownames(A5), rownames(A6), rownames(A7), rownames(A8), rownames(A9), rownames(A10)))

### only merge samples present in all
A3<-A3[probes_all,]
A4<-A4[probes_all,]
A5<-A5[probes_all,]
A6<-A6[probes_all,]
A7<-A7[probes_all,]
A8<-A8[probes_all,]
A9<-A9[probes_all,]
A10<-A10[probes_all,]

### summary statistics per sample:

summary.id<-matrix(data = NA, nrow = 8, ncol = 11)
colnames(summary.id)<-c("Number of Sites", "Min Read Depth", "1st Quantile Read Depth", "Mean Read Depth", "Median Read Depth",  "3rd Quantile Read Depth", "Max Read Depth", "Percentage Sites:Read Depth 5",  "Percentage Sites:Read Depth 10", "Sites:Read Depth 5","Sites:Read Depth 10")
rownames(summary.id)<-c("A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")

summary.id[1,1]<-nrow(A3)
summary.id[2,1]<-nrow(A4)
summary.id[3,1]<-nrow(A5)
summary.id[4,1]<-nrow(A6)
summary.id[5,1]<-nrow(A7)
summary.id[6,1]<-nrow(A8)
summary.id[7,1]<-nrow(A9)
summary.id[8,1]<-nrow(A10)

summary.id[1,2]<-min(A3[,4])
summary.id[2,2]<-min(A4[,4])
summary.id[3,2]<-min(A5[,4])
summary.id[4,2]<-min(A6[,4])
summary.id[5,2]<-min(A7[,4])
summary.id[6,2]<-min(A8[,4])
summary.id[7,2]<-min(A9[,4])
summary.id[8,2]<-min(A10[,4])

summary.id[1,3]<-quantile(A3[,4], 0.25)
summary.id[2,3]<-quantile(A4[,4], 0.25)
summary.id[3,3]<-quantile(A5[,4], 0.25)
summary.id[4,3]<-quantile(A6[,4], 0.25)
summary.id[5,3]<-quantile(A7[,4], 0.25)
summary.id[6,3]<-quantile(A8[,4], 0.25)
summary.id[7,3]<-quantile(A9[,4], 0.25)
summary.id[8,3]<-quantile(A10[,4], 0.25)

summary.id[1,4]<-mean(A3[,4])
summary.id[2,4]<-mean(A4[,4])
summary.id[3,4]<-mean(A5[,4])
summary.id[4,4]<-mean(A6[,4])
summary.id[5,4]<-mean(A7[,4])
summary.id[6,4]<-mean(A8[,4])
summary.id[7,4]<-mean(A9[,4])
summary.id[8,4]<-mean(A10[,4])

summary.id[1,5]<-median(A3[,4])
summary.id[2,5]<-median(A4[,4])
summary.id[3,5]<-median(A5[,4])
summary.id[4,5]<-median(A6[,4])
summary.id[5,5]<-median(A7[,4])
summary.id[6,5]<-median(A8[,4])
summary.id[7,5]<-median(A9[,4])
summary.id[8,5]<-median(A10[,4])

summary.id[1,6]<-quantile(A3[,4], 0.75)
summary.id[2,6]<-quantile(A4[,4], 0.75)
summary.id[3,6]<-quantile(A5[,4], 0.75)
summary.id[4,6]<-quantile(A6[,4], 0.75)
summary.id[5,6]<-quantile(A7[,4], 0.75)
summary.id[6,6]<-quantile(A8[,4], 0.75)
summary.id[7,6]<-quantile(A9[,4], 0.75)
summary.id[8,6]<-quantile(A10[,4], 0.75)

summary.id[1,7]<-max(A3[,4])
summary.id[2,7]<-max(A4[,4])
summary.id[3,7]<-max(A5[,4])
summary.id[4,7]<-max(A6[,4])
summary.id[5,7]<-max(A7[,4])
summary.id[6,7]<-max(A8[,4])
summary.id[7,7]<-max(A9[,4])
summary.id[8,7]<-max(A10[,4])

summary.id[1,8]<-length(which(A3[,4] >= 5))/nrow(A3)*100
summary.id[2,8]<-length(which(A4[,4] >= 5))/nrow(A4)*100
summary.id[3,8]<-length(which(A5[,4] >= 5))/nrow(A5)*100
summary.id[4,8]<-length(which(A6[,4] >= 5))/nrow(A6)*100
summary.id[5,8]<-length(which(A7[,4] >= 5))/nrow(A7)*100
summary.id[6,8]<-length(which(A8[,4] >= 5))/nrow(A8)*100
summary.id[7,8]<-length(which(A9[,4] >= 5))/nrow(A9)*100
summary.id[8,8]<-length(which(A10[,4] >= 5))/nrow(A10)*100

summary.id[1,9]<-length(which(A3[,4] >= 10))/nrow(A3)*100
summary.id[2,9]<-length(which(A4[,4] >= 10))/nrow(A4)*100
summary.id[3,9]<-length(which(A5[,4] >= 10))/nrow(A5)*100
summary.id[4,9]<-length(which(A6[,4] >= 10))/nrow(A6)*100
summary.id[5,9]<-length(which(A7[,4] >= 10))/nrow(A7)*100
summary.id[6,9]<-length(which(A8[,4] >= 10))/nrow(A8)*100
summary.id[7,9]<-length(which(A9[,4] >= 10))/nrow(A9)*100
summary.id[8,9]<-length(which(A10[,4] >= 10))/nrow(A10)*100

summary.id[1,10]<-length(which(A3[,4] >= 5))
summary.id[2,10]<-length(which(A4[,4] >= 5))
summary.id[3,10]<-length(which(A5[,4] >= 5))
summary.id[4,10]<-length(which(A6[,4] >= 5))
summary.id[5,10]<-length(which(A7[,4] >= 5))
summary.id[6,10]<-length(which(A8[,4] >= 5))
summary.id[7,10]<-length(which(A9[,4] >= 5))
summary.id[8,10]<-length(which(A10[,4] >= 5))

summary.id[1,11]<-length(which(A3[,4] >= 10))
summary.id[2,11]<-length(which(A4[,4] >= 10))
summary.id[3,11]<-length(which(A5[,4] >= 10))
summary.id[4,11]<-length(which(A6[,4] >= 10))
summary.id[5,11]<-length(which(A7[,4] >= 10))
summary.id[6,11]<-length(which(A8[,4] >= 10))
summary.id[7,11]<-length(which(A9[,4] >= 10))
summary.id[8,11]<-length(which(A10[,4] >= 10))

write.csv(summary.id, "EachSampleReadDepthStatsCommonSites_CpGsites.csv")

all<-cbind(A3,A4,A5,A6,A7,A8,A9,A10)

## calc min, median, mean and max read depth across samples

coveragecolumns<-grep("TotReads", colnames(all))

CalcnSamplesPres<-function(i){
	return(length(which(is.na(i) == FALSE)))

}

nPres<-apply(all[,coveragecolumns], 1, CalcnSamplesPres)
minCov<-apply(all[,coveragecolumns], 1, min, na.rm = T)
medCov<-apply(all[,coveragecolumns], 1, median, na.rm = T)
meanCov<-apply(all[,coveragecolumns], 1, mean, na.rm = T)
maxCov<-apply(all[,coveragecolumns], 1, max, na.rm = T)

all<-cbind(all, nPres, minCov,medCov,meanCov,maxCov)
### Filter to sites present in all:

write.csv(all, "MethylationData_CpG_PresentAllSamples.csv")
save(all, file = "CpG_PresentAllSamples.Rdata")


### read statistics per across sites

summary.site<-matrix(data = NA, nrow = 4, ncol = 12)
rownames(summary.site)<-c("Site in any sample", "Min Read Depth 1", "Min Read Depth 5", "Min Read Depth 10") 
colnames(summary.site)<-c("Number of Sites", "Percentage of All", "Min Min Read Depth", "1st Quantile Min Read Depth", "Mean Min Read Depth", "Median Min Read Depth",  "3rd Quantile Min Read Depth", "Max Min Read Depth", "Percentage sites with Min Read Depth 5", "Percentage sites with Min Read Depth 10","Percentage sites with Min Read Depth 15","Percentage sites with Min Read Depth 20")

summary.site[1,1]<-length(all_poss)

summary.site[2,1]<-nrow(all)
summary.site[2,2]<-nrow(all)/nrow(all)*100
summary.site[2,3]<-min(all[,34])
summary.site[2,4]<-quantile(all[,34], 0.25)
summary.site[2,5]<-mean(all[,34])
summary.site[2,6]<-median(all[,34])
summary.site[2,7]<-quantile(all[,34], 0.75)
summary.site[2,8]<-max(all[,34])
summary.site[2,9]<-length(which(all[,34] >= 5))/nrow(all)*100
summary.site[2,10]<-length(which(all[,34] >= 10))/nrow(all)*100
summary.site[2,11]<-length(which(all[,34] >= 15))/nrow(all)*100
summary.site[2,12]<-length(which(all[,34] >= 20))/nrow(all)*100

all.sub<-subset(all, minCov >= 5)
summary.site[3,1]<-nrow(all.sub)
summary.site[3,2]<-nrow(all.sub)/nrow(all)*100
summary.site[3,3]<-min(all.sub[,34])
summary.site[3,4]<-quantile(all.sub[,34], 0.25)
summary.site[3,5]<-mean(all.sub[,34])
summary.site[3,6]<-median(all.sub[,34])
summary.site[3,7]<-quantile(all.sub[,34], 0.75)
summary.site[3,8]<-max(all.sub[,34])
summary.site[3,9]<-length(which(all.sub[,34] >= 5))/nrow(all.sub)*100
summary.site[3,10]<-length(which(all.sub[,34] >= 10))/nrow(all.sub)*100
summary.site[3,11]<-length(which(all.sub[,34] >= 15))/nrow(all.sub)*100
summary.site[3,12]<-length(which(all.sub[,34] >= 20))/nrow(all.sub)*100

all.sub<-subset(all, minCov >= 10)
summary.site[4,1]<-nrow(all.sub)
summary.site[4,2]<-nrow(all.sub)/nrow(all)*100
summary.site[4,3]<-min(all.sub[,34])
summary.site[4,4]<-quantile(all.sub[,34], 0.25)
summary.site[4,5]<-mean(all.sub[,34])
summary.site[4,6]<-median(all.sub[,34])
summary.site[4,7]<-quantile(all.sub[,34], 0.75)
summary.site[4,8]<-max(all.sub[,34])
summary.site[4,9]<-length(which(all.sub[,34] >= 5))/nrow(all.sub)*100
summary.site[4,10]<-length(which(all.sub[,34] >= 10))/nrow(all.sub)*100
summary.site[4,11]<-length(which(all.sub[,34] >= 15))/nrow(all.sub)*100
summary.site[4,12]<-length(which(all.sub[,34] >= 20))/nrow(all.sub)*100

write.csv(summary.site, "SiteReadDepthStatsCommonSites_CpGsites.csv")

### boxplot of read depth

pdf("/mnt/data1/Eilis/Projects/Rat_Neuron/plots/BoxplotReadDepth_CpG_sitesPresentinallSamples.pdf")
boxplot(list(all[,4], all[,8],all[,12],all[,16],all[,20],all[,24],all[,28],all[,32]), ylab = "Read Depth", main = "All sites present in all Samples", names = c("A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"))
boxplot(list(all[,4], all[,8],all[,12],all[,16],all[,20],all[,24],all[,28],all[,32]), ylab = "Read Depth", main = "All sites present in all Samples", names = c("A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"), ylim = c(0,1000))
boxplot(list(all[,34], all[,36], all[,35], all[,37]), ylab = "Read Depth", xlab = "Read Depth across all Samples", main = "All sites present in all Samples", names = c("Minimum", "Mean", "Median", "Maximum"))
boxplot(list(all[,34], all[,36], all[,35], all[,37]), ylab = "Read Depth", xlab = "Read Depth across all Samples", main = "All sites present in all Samples", names = c("Minimum", "Mean", "Median", "Maximum"), ylim = c(0,1000))
dev.off()
### Identify most variable probes
methcols<-grep("Methylation", colnames(all))

methsigma<-apply(all[, methcols],1, sd)
methmu<-apply(all[, methcols],1, mean)

tab<-matrix(data = NA, nrow = 8, ncol = 8)
for(i in 1:8){
	for(j in 1:8){
		tab[i,j]<-cor(all[,methcols[i]], all[, methcols[j]], method = "pearson")
	}
}

sigmarank<-rank(methsigma)
topVariableProbes<-which(sigmarank > (length(methsigma)-10000))
d<-dist(t(all[topVariableProbes, methcols]))
h<-hclust(d)
pdf("HierarchicalClusterofSamples2_CpG.pdf")
plot(h, main = "Top 10001 Variable Probes", labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
fit <- cmdscale(d,eig=TRUE, k=2)

plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c("blue", "red", "blue", "red", "blue", "red", "blue", "red")) 
legend("topright", c("Control", "K"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c(rep("blue",4), rep("red",4))) 
legend("topright", c("June11", "Oct11"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "MDS Top 10001 Variable Probes", col = c("blue", "blue", "red", "red", "blue", "blue", "green", "green"), pch = c(2,3,2,3,2,3,2,3))
legend("topright", c("NA", "nifedipine", "FK-506", "Control", "K"), col = c("blue", "red", "green", "black", "black"), pch = c(1,1,1,2,3))
boxplot(all[, methcols], las = 2, labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
dev.off()


### Cluster on most variable sites with min read depth of 10
all.10<-subset(all, minCov > 9)
methsigma<-apply(all.10[, methcols],1, sd)
sigmarank<-rank(methsigma)
topVariableProbes<-which(sigmarank > (length(methsigma)-10000))
d<-dist(t(all.10[topVariableProbes, methcols]))
h<-hclust(d)
pdf("HierarchicalClusterofSamplesminReadDepth10_CpG.pdf")
plot(h, main = "Top 10001 Variable Probes with min Read Depth 10", labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))

### MDS of same probes
fit <- cmdscale(d,eig=TRUE, k=2)

plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c("blue", "red", "blue", "red", "blue", "red", "blue", "red")) 
legend("topright", c("Control", "K"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c(rep("blue",4), rep("red",4))) 
legend("topright", c("June11", "Oct11"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "MDS Top 10001 Variable Probes", col = c("blue", "blue", "red", "red", "blue", "blue", "green", "green"), pch = c(2,3,2,3,2,3,2,3))
legend("topright", c("NA", "nifedipine", "FK-506", "Control", "K"), col = c("blue", "red", "green", "black", "black"), pch = c(1,1,1,2,3))
boxplot(all.10[, methcols], las = 2, labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
dev.off()


### exclude sex chromosomes
all<-all[-grep("chrX", rownames(all)),]
all<-all[-grep("chrM", rownames(all)),]

methcols<-grep("Methylation", colnames(all))
methsigma<-apply(all[, methcols],1, sd)
sigmarank<-rank(methsigma)
topVariableProbes<-which(sigmarank > (length(methsigma)-10000))
d<-dist(t(all[topVariableProbes, methcols]))
h<-hclust(d)
pdf("HierarchicalClusterofSamples2_CpG_excludeXChr.pdf")
plot(h, main = "Top 10001 Variable Probes", labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
fit <- cmdscale(d,eig=TRUE, k=2)

plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c("blue", "red", "blue", "red", "blue", "red", "blue", "red")) 
legend("topleft", c("Control", "K"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c(rep("blue",4), rep("red",4))) 
legend("topleft", c("June11", "Oct11"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "MDS Top 10001 Variable Probes", col = c("blue", "blue", "red", "red", "blue", "blue", "green", "green"), pch = c(2,3,2,3,2,3,2,3))
legend("topleft", c("NA", "nifedipine", "FK-506", "Control", "K"), col = c("blue", "red", "green", "black", "black"), pch = c(1,1,1,2,3))
boxplot(all[, methcols], las = 2, labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
dev.off()


### Cluster on most variable sites with min read depth of 10
all.10<-subset(all, minCov > 9)
methsigma<-apply(all.10[, methcols],1, sd)
sigmarank<-rank(methsigma)
topVariableProbes<-which(sigmarank > (length(methsigma)-10000))
d<-dist(t(all.10[topVariableProbes, methcols]))
h<-hclust(d)
pdf("HierarchicalClusterofSamplesminReadDepth10_CpG_ExcludeXchr.pdf")
plot(h, main = "Top 10001 Variable Probes with min Read Depth 10", labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))

### MDS of same probes
fit <- cmdscale(d,eig=TRUE, k=2)

plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c("blue", "red", "blue", "red", "blue", "red", "blue", "red")) 
legend("topleft", c("Control", "K"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", type = "n", main = "MDS Top 10001 Variable Probes")
text(fit$points[,1], fit$points[,2], labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"), col = c(rep("blue",4), rep("red",4))) 
legend("topleft", c("June11", "Oct11"), col = c("blue", "red"), pch = 1)
plot(fit$points[,1], fit$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", main = "MDS Top 10001 Variable Probes", col = c("blue", "blue", "red", "red", "blue", "blue", "green", "green"), pch = c(2,3,2,3,2,3,2,3))
legend("topleft", c("NA", "nifedipine", "FK-506", "Control", "K"), col = c("blue", "red", "green", "black", "black"), pch = c(1,1,1,2,3))
boxplot(all.10[, methcols], las = 2, labels = c("A3","A4", "A5","A6","A7", "A8", "A9", "A10"))
dev.off()

