#!/bin/bash

### To run script:
# add permissions to execute: chmod +x 02-mapping_RRBS_pilot.sh
# sh 02-mapping_RRBS_pilot.sh

# Isabel Castanho # I.Castanho@exeter.ac.uk
# RRBS
# alignment (3-letter alignment of bisulfite-Seq reads) using Bismark (https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)


# FastQC
# done by the Seq service


# reference genome: /mnt/data1/reference_files/mouse_GRCm38.p5/GCA_000001635.7_GRCm38.p5_genomic.fa
# spike-in sequences appended to the ref genome fasta file by Joana: GCA_000001635.7_GRCm38.p5_genomic.spikein.fa


# Trimming using Trim galore
# remember that the adapters are added before bisulfite conversion, so any C's need to be T's now (from Joana)
# change error rate to 0.2 to allow at least 1 mismatch
cd /mnt/data1/isabel/RRBS/01_miseq_nano/raw_reads
input="/mnt/data1/isabel/RRBS/01_miseq_nano/raw_reads/*_r1.fq.gz" # If the input files are gzip-compressed the output files will be automatically gzip compressed as well.
output="/mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/"

# /mnt/data1/programs/trim_galore_v0.3.3/trim_galore # old version! Update next time!
perl /home/icastanho/trim_galore --quality 20 --phred33 --rrbs $input --output_dir $output --fastqc_args "--outdir /mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/FASTQC/" 
# use --ilumina in future versions
# I can now run just typing trim_galore in the command line (16-Dec-2017)
# cutadapt can know be called from the command line as well


# (I) Running bismark_genome_preparation # run bismark on ref genome + spike-in # directional???
/mnt/data1/programs/Bismark/bismark_genome_preparation --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1 --verbose /mnt/data1/isabel/RRBS/01_miseq_nano/reference_genome/

# (II) Running bismark # 
for file in ./*_r1_trimmed.fq.gz
do
	echo $file
	bismark -N 1 --un --ambiguous --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1/ --genome /mnt/data1/isabel/RRBS/01_miseq_nano/reference_genome/ --se $file
done

# convert .bam to .sam
#for file in ./*.bam
#do
#	echo $file 
#	samtools view -h $file > ${file/.bam/.sam}
#done

# sort and index files
# sort
#for file in ./*.bam
#do
#	echo $file
#	samtools sort $file > $file.sorted.bam
#done

# index
#for file in ./*.sorted.bam
#do
#	echo $file
#	samtools index $file > $file.sorted.indexed.bam
#done


# SeqMonk
# /mnt/data1/programs/SeqMonk/seqmonk # v 1.39.0 (updated 19-Dec-2017)
/mnt/data1/programs/SeqMonk/seqmonk

# (III) Running bismark_methylation_extractor 
for file in ./*_r1_trimmed_bismark_bt2.bam
do
	echo $file
	bismark_methylation_extractor -s --gzip --comprehensive --multicore 4 --bedGraph --cytosine_report --genome_folder /mnt/data1/isabel/RRBS/01_miseq_nano/reference_genome/ $file
done


# (IV) Running bismark2report
bismark2report

# (V) Running bismark2summary
bismark2summary # bismark_summary_report.html

# convert BAM to SAM
for file in ./*.bam
	do
	echo $file 
	samtools view -h $file > ${file/.bam/.sam}
	done

### create unique lists of methylation sites split by type # by Eilis
zcat CpG_context_*_r1_trimmed_bismark_bt2.txt.gz | cut -f3,4 | sort | uniq -c > All_CpG_Sites_AcrossAllSamples.txt
zcat CHG_context_*_r1_trimmed_bismark_bt2.txt.gz | cut -f3,4 | sort | uniq -c > All_CHG_Sites_AcrossAllSamples.txt
zcat CHH_context_*_r1_trimmed_bismark_bt2.txt.gz | cut -f3,4 | sort | uniq -c > All_CHH_Sites_AcrossAllSamples.txt

### create list of all sites found in at least one sample
zcat *_r1_trimmed_bismark_bt2.bismark.cov.gz | cut -f 1-2 | sort | uniq -c | wc


### calculate sumary statistics separately for CpG and nonCpG methylation
R

setwd("/mnt/data1/isabel/RRBS/01_miseq_nano/methylation_extraction/")

#### From Eilis
### start with non CpG
cpg <- read.table("All_CpG_Sites_AcrossAllSamples.txt", header=TRUE, fill=TRUE)
cpg <- paste(cpg[,2], cpg[,3], sep = "_")

##create list of probes for each file

CreateProbeList<-function(x){
  rownames(x)<-paste(x[,1], x[,2], sep = "_")
  return(x)
}

##Load Bedgraph files
A17<-read.table("A17_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A17<-CreateProbeList(A17)
A17<-A17[intersect(cpg, rownames(A17)),]
A18<-read.table("A24_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A18<-CreateProbeList(A18)
A18<-A18[intersect(cpg, rownames(A18)),]
A19<-read.table("A26_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A19<-CreateProbeList(A19)
A19<-A19[intersect(cpg, rownames(A19)),]
A20<-read.table("A28_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A20<-CreateProbeList(A20)
A20<-A20[intersect(cpg, rownames(A20)),]
A21<-read.table("B29_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A21<-CreateProbeList(A21)
A21<-A21[intersect(cpg, rownames(A21)),]
A22<-read.table("B30_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A22<-CreateProbeList(A22)
A22<-A22[intersect(cpg, rownames(A22)),]
A23<-read.table("B31_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A23<-CreateProbeList(A23)
A23<-A23[intersect(cpg, rownames(A23)),]
A24<-read.table("B32_r1_trimmed_bismark_bt2.bismark.cov", header=FALSE)
A24<-CreateProbeList(A24)
A24<-A24[intersect(cpg, rownames(A24)),]

A17<-as.matrix(A17[,-c(1:3)])
A18<-as.matrix(A18[,-c(1:3)])
A19<-as.matrix(A19[,-c(1:3)])
A20<-as.matrix(A20[,-c(1:3)])
A21<-as.matrix(A21[,-c(1:3)])
A22<-as.matrix(A22[,-c(1:3)])
A23<-as.matrix(A23[,-c(1:3)])
A24<-as.matrix(A24[,-c(1:3)])

A17<-A17[intersect(cpg, rownames(A17)),]
A18<-A18[intersect(cpg, rownames(A18)),]
A19<-A19[intersect(cpg, rownames(A19)),]
A20<-A20[intersect(cpg, rownames(A20)),]
A21<-A21[intersect(cpg, rownames(A21)),]
A22<-A22[intersect(cpg, rownames(A22)),]
A23<-A23[intersect(cpg, rownames(A23)),]
A24<-A24[intersect(cpg, rownames(A24)),]

save(A17,A18,A19,A20,A21,A22,A23,A24, file="CpGMethylation.Rdata")

### calculate coverage:
CalcCoverage<-function(x, lab){
  x<-cbind(x, x[,2]+x[,3])
  colnames(x)<-paste(lab, c('Methylation', 'nMReads', 'nNMReads', 'TotReads'), sep = "_")
  return(x)
}

A17<-CalcCoverage(A17, lab = "A17")
A18<-CalcCoverage(A18, lab = "A18")
A19<-CalcCoverage(A19, lab = "A19")
A20<-CalcCoverage(A20, lab = "A20")
A21<-CalcCoverage(A21, lab = "A21")
A22<-CalcCoverage(A22, lab = "A22")
A23<-CalcCoverage(A23, lab = "A23")
A24<-CalcCoverage(A24, lab = "A24")

save(A17,A18,A19,A20,A21,A22,A23,A24, file = "CpGMethylation.Rdata")

### summary statistics per sample:

summary.id<-matrix(data = NA, nrow = 8, ncol = 11)
colnames(summary.id)<-c("Number of Sites", "Min Read Depth", "1st Quantile Read Depth", "Mean Read Depth", "Median Read Depth",  "3rd Quantile Read Depth", "Max Read Depth", "Percentage Sites:Read Depth 5",  "Percentage Sites:Read Depth 10", "Sites:Read Depth 5","Sites:Read Depth 10")
rownames(summary.id)<-c("A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24")

summary.id[1,1]<-nrow(A17)
summary.id[2,1]<-nrow(A18)
summary.id[3,1]<-nrow(A19)
summary.id[4,1]<-nrow(A20)
summary.id[5,1]<-nrow(A21)
summary.id[6,1]<-nrow(A22)
summary.id[7,1]<-nrow(A23)
summary.id[8,1]<-nrow(A24)

summary.id[1,2]<-min(A17[,4])
summary.id[2,2]<-min(A18[,4])
summary.id[3,2]<-min(A19[,4])
summary.id[4,2]<-min(A20[,4])
summary.id[5,2]<-min(A21[,4])
summary.id[6,2]<-min(A22[,4])
summary.id[7,2]<-min(A23[,4])
summary.id[8,2]<-min(A24[,4])

summary.id[1,3]<-quantile(A17[,4], 0.25)
summary.id[2,3]<-quantile(A18[,4], 0.25)
summary.id[3,3]<-quantile(A19[,4], 0.25)
summary.id[4,3]<-quantile(A20[,4], 0.25)
summary.id[5,3]<-quantile(A21[,4], 0.25)
summary.id[6,3]<-quantile(A22[,4], 0.25)
summary.id[7,3]<-quantile(A23[,4], 0.25)
summary.id[8,3]<-quantile(A24[,4], 0.25)

summary.id[1,4]<-mean(A17[,4])
summary.id[2,4]<-mean(A18[,4])
summary.id[3,4]<-mean(A19[,4])
summary.id[4,4]<-mean(A20[,4])
summary.id[5,4]<-mean(A21[,4])
summary.id[6,4]<-mean(A22[,4])
summary.id[7,4]<-mean(A23[,4])
summary.id[8,4]<-mean(A24[,4])

summary.id[1,5]<-median(A17[,4])
summary.id[2,5]<-median(A18[,4])
summary.id[3,5]<-median(A19[,4])
summary.id[4,5]<-median(A20[,4])
summary.id[5,5]<-median(A21[,4])
summary.id[6,5]<-median(A22[,4])
summary.id[7,5]<-median(A23[,4])
summary.id[8,5]<-median(A24[,4])

summary.id[1,6]<-quantile(A17[,4], 0.75)
summary.id[2,6]<-quantile(A18[,4], 0.75)
summary.id[3,6]<-quantile(A19[,4], 0.75)
summary.id[4,6]<-quantile(A20[,4], 0.75)
summary.id[5,6]<-quantile(A21[,4], 0.75)
summary.id[6,6]<-quantile(A22[,4], 0.75)
summary.id[7,6]<-quantile(A23[,4], 0.75)
summary.id[8,6]<-quantile(A24[,4], 0.75)

summary.id[1,7]<-max(A17[,4])
summary.id[2,7]<-max(A18[,4])
summary.id[3,7]<-max(A19[,4])
summary.id[4,7]<-max(A20[,4])
summary.id[5,7]<-max(A21[,4])
summary.id[6,7]<-max(A22[,4])
summary.id[7,7]<-max(A23[,4])
summary.id[8,7]<-max(A24[,4])

summary.id[1,8]<-length(which(A17[,4] >= 5))/nrow(A17)*100
summary.id[2,8]<-length(which(A18[,4] >= 5))/nrow(A18)*100
summary.id[3,8]<-length(which(A19[,4] >= 5))/nrow(A19)*100
summary.id[4,8]<-length(which(A20[,4] >= 5))/nrow(A20)*100
summary.id[5,8]<-length(which(A21[,4] >= 5))/nrow(A21)*100
summary.id[6,8]<-length(which(A22[,4] >= 5))/nrow(A22)*100
summary.id[7,8]<-length(which(A23[,4] >= 5))/nrow(A23)*100
summary.id[8,8]<-length(which(A24[,4] >= 5))/nrow(A24)*100

summary.id[1,9]<-length(which(A17[,4] >= 10))/nrow(A17)*100
summary.id[2,9]<-length(which(A18[,4] >= 10))/nrow(A18)*100
summary.id[3,9]<-length(which(A19[,4] >= 10))/nrow(A19)*100
summary.id[4,9]<-length(which(A20[,4] >= 10))/nrow(A20)*100
summary.id[5,9]<-length(which(A21[,4] >= 10))/nrow(A21)*100
summary.id[6,9]<-length(which(A22[,4] >= 10))/nrow(A22)*100
summary.id[7,9]<-length(which(A23[,4] >= 10))/nrow(A23)*100
summary.id[8,9]<-length(which(A24[,4] >= 10))/nrow(A24)*100

summary.id[1,10]<-length(which(A17[,4] >= 5))
summary.id[2,10]<-length(which(A18[,4] >= 5))
summary.id[3,10]<-length(which(A19[,4] >= 5))
summary.id[4,10]<-length(which(A20[,4] >= 5))
summary.id[5,10]<-length(which(A21[,4] >= 5))
summary.id[6,10]<-length(which(A22[,4] >= 5))
summary.id[7,10]<-length(which(A23[,4] >= 5))
summary.id[8,10]<-length(which(A24[,4] >= 5))

summary.id[1,11]<-length(which(A17[,4] >= 10))
summary.id[2,11]<-length(which(A18[,4] >= 10))
summary.id[3,11]<-length(which(A19[,4] >= 10))
summary.id[4,11]<-length(which(A20[,4] >= 10))
summary.id[5,11]<-length(which(A21[,4] >= 10))
summary.id[6,11]<-length(which(A22[,4] >= 10))
summary.id[7,11]<-length(which(A23[,4] >= 10))
summary.id[8,11]<-length(which(A24[,4] >= 10))

write.csv(summary.id, "EachSampleReadDepthStatsAllSites_CpGsites.csv")

### Create List of probes covered in all samples
probes_all<-intersect(rownames(A17), rownames(A18))
probes_all<-intersect(probes_all, rownames(A19))
probes_all<-intersect(probes_all, rownames(A20))
probes_all<-intersect(probes_all, rownames(A21))
probes_all<-intersect(probes_all, rownames(A22))
probes_all<-intersect(probes_all, rownames(A23))
probes_all<-intersect(probes_all, rownames(A24))

##obtain total number of sites assessed across all samples
allposs<-unique(c(rownames(A17), rownames(A18), rownames(A19), rownames(A20), rownames(A21), rownames(A22), rownames(A23), rownames(A24)))

### only merge samples present in all
A17<-A17[probes_all,]
A18<-A18[probes_all,]
A19<-A19[probes_all,]
A20<-A20[probes_all,]
A21<-A21[probes_all,]
A22<-A22[probes_all,]
A23<-A23[probes_all,]
A24<-A24[probes_all,]

### summary statistics per sample:

summary.id<-matrix(data = NA, nrow = 8, ncol = 11)
colnames(summary.id)<-c("Number of Sites", "Min Read Depth", "1st Quantile Read Depth", "Mean Read Depth", "Median Read Depth",  "3rd Quantile Read Depth", "Max Read Depth", "Percentage Sites:Read Depth 5",  "Percentage Sites:Read Depth 10", "Sites:Read Depth 5","Sites:Read Depth 10")
rownames(summary.id)<-c("A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24")

summary.id[1,1]<-nrow(A17)
summary.id[2,1]<-nrow(A18)
summary.id[3,1]<-nrow(A19)
summary.id[4,1]<-nrow(A20)
summary.id[5,1]<-nrow(A21)
summary.id[6,1]<-nrow(A22)
summary.id[7,1]<-nrow(A23)
summary.id[8,1]<-nrow(A24)

summary.id[1,2]<-min(A17[,4])
summary.id[2,2]<-min(A18[,4])
summary.id[3,2]<-min(A19[,4])
summary.id[4,2]<-min(A20[,4])
summary.id[5,2]<-min(A21[,4])
summary.id[6,2]<-min(A22[,4])
summary.id[7,2]<-min(A23[,4])
summary.id[8,2]<-min(A24[,4])

summary.id[1,3]<-quantile(A17[,4], 0.25)
summary.id[2,3]<-quantile(A18[,4], 0.25)
summary.id[3,3]<-quantile(A19[,4], 0.25)
summary.id[4,3]<-quantile(A20[,4], 0.25)
summary.id[5,3]<-quantile(A21[,4], 0.25)
summary.id[6,3]<-quantile(A22[,4], 0.25)
summary.id[7,3]<-quantile(A23[,4], 0.25)
summary.id[8,3]<-quantile(A24[,4], 0.25)

summary.id[1,4]<-mean(A17[,4])
summary.id[2,4]<-mean(A18[,4])
summary.id[3,4]<-mean(A19[,4])
summary.id[4,4]<-mean(A20[,4])
summary.id[5,4]<-mean(A21[,4])
summary.id[6,4]<-mean(A22[,4])
summary.id[7,4]<-mean(A23[,4])
summary.id[8,4]<-mean(A24[,4])

summary.id[1,5]<-median(A17[,4])
summary.id[2,5]<-median(A18[,4])
summary.id[3,5]<-median(A19[,4])
summary.id[4,5]<-median(A20[,4])
summary.id[5,5]<-median(A21[,4])
summary.id[6,5]<-median(A22[,4])
summary.id[7,5]<-median(A23[,4])
summary.id[8,5]<-median(A24[,4])

summary.id[1,6]<-quantile(A17[,4], 0.75)
summary.id[2,6]<-quantile(A18[,4], 0.75)
summary.id[3,6]<-quantile(A19[,4], 0.75)
summary.id[4,6]<-quantile(A20[,4], 0.75)
summary.id[5,6]<-quantile(A21[,4], 0.75)
summary.id[6,6]<-quantile(A22[,4], 0.75)
summary.id[7,6]<-quantile(A23[,4], 0.75)
summary.id[8,6]<-quantile(A24[,4], 0.75)

summary.id[1,7]<-max(A17[,4])
summary.id[2,7]<-max(A18[,4])
summary.id[3,7]<-max(A19[,4])
summary.id[4,7]<-max(A20[,4])
summary.id[5,7]<-max(A21[,4])
summary.id[6,7]<-max(A22[,4])
summary.id[7,7]<-max(A23[,4])
summary.id[8,7]<-max(A24[,4])

summary.id[1,8]<-length(which(A17[,4] >= 5))/nrow(A17)*100
summary.id[2,8]<-length(which(A18[,4] >= 5))/nrow(A18)*100
summary.id[3,8]<-length(which(A19[,4] >= 5))/nrow(A19)*100
summary.id[4,8]<-length(which(A20[,4] >= 5))/nrow(A20)*100
summary.id[5,8]<-length(which(A21[,4] >= 5))/nrow(A21)*100
summary.id[6,8]<-length(which(A22[,4] >= 5))/nrow(A22)*100
summary.id[7,8]<-length(which(A23[,4] >= 5))/nrow(A23)*100
summary.id[8,8]<-length(which(A24[,4] >= 5))/nrow(A24)*100

summary.id[1,9]<-length(which(A17[,4] >= 10))/nrow(A17)*100
summary.id[2,9]<-length(which(A18[,4] >= 10))/nrow(A18)*100
summary.id[3,9]<-length(which(A19[,4] >= 10))/nrow(A19)*100
summary.id[4,9]<-length(which(A20[,4] >= 10))/nrow(A20)*100
summary.id[5,9]<-length(which(A21[,4] >= 10))/nrow(A21)*100
summary.id[6,9]<-length(which(A22[,4] >= 10))/nrow(A22)*100
summary.id[7,9]<-length(which(A23[,4] >= 10))/nrow(A23)*100
summary.id[8,9]<-length(which(A24[,4] >= 10))/nrow(A24)*100

summary.id[1,10]<-length(which(A17[,4] >= 5))
summary.id[2,10]<-length(which(A18[,4] >= 5))
summary.id[3,10]<-length(which(A19[,4] >= 5))
summary.id[4,10]<-length(which(A20[,4] >= 5))
summary.id[5,10]<-length(which(A21[,4] >= 5))
summary.id[6,10]<-length(which(A22[,4] >= 5))
summary.id[7,10]<-length(which(A23[,4] >= 5))
summary.id[8,10]<-length(which(A24[,4] >= 5))

summary.id[1,11]<-length(which(A17[,4] >= 10))
summary.id[2,11]<-length(which(A18[,4] >= 10))
summary.id[3,11]<-length(which(A19[,4] >= 10))
summary.id[4,11]<-length(which(A20[,4] >= 10))
summary.id[5,11]<-length(which(A21[,4] >= 10))
summary.id[6,11]<-length(which(A22[,4] >= 10))
summary.id[7,11]<-length(which(A23[,4] >= 10))
summary.id[8,11]<-length(which(A24[,4] >= 10))

write.csv(summary.id, "EachSampleReadDepthStatsCommonSites_CpGsites.csv")

all<-cbind(A17,A18,A19,A20,A21,A22,A23,A24)

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
pdf("/mnt/data1/isabel/RRBS/01_miseq_nano/methylation_extraction/BoxplotReadDepth_CpG_sitesPresentinallSamples.pdf")
boxplot(list(all[,4], all[,8],all[,12],all[,16],all[,20],all[,24],all[,28],all[,32]), ylab = "Read Depth", main = "All sites present in all Samples", names = c("A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24"))
boxplot(list(all[,4], all[,8],all[,12],all[,16],all[,20],all[,24],all[,28],all[,32]), ylab = "Read Depth", main = "All sites present in all Samples", names = c("A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24"), ylim = c(-10,300))
boxplot(list(all[,4], all[,8],all[,12],all[,16],all[,20],all[,24],all[,28],all[,32]), ylab = "Read Depth", main = "All sites present in all Samples", names = c("A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24"), ylim = c(-10,100))
boxplot(list(all[,34], all[,36], all[,35], all[,37]), ylab = "Read Depth", xlab = "Read Depth across all Samples", main = "All sites present in all Samples", names = c("Minimum", "Mean", "Median", "Maximum"))
boxplot(list(all[,34], all[,36], all[,35], all[,37]), ylab = "Read Depth", xlab = "Read Depth across all Samples", main = "All sites present in all Samples", names = c("Minimum", "Mean", "Median", "Maximum"), ylim = c(-10,300))
boxplot(list(all[,34], all[,36], all[,35], all[,37]), ylab = "Read Depth", xlab = "Read Depth across all Samples", main = "All sites present in all Samples", names = c("Minimum", "Mean", "Median", "Maximum"), ylim = c(-10,100))
dev.off()

### Histograms of beta values for NonCpg methylatio (together and sep as CHG and CHH)
load("CpG_PresentAllSamples.Rdata")

### filter to min read depth of 10
sub_all<-all[which(all[,34] >= 10),]

### histogram of beta values per sample
pdf("/mnt/data1/isabel/RRBS/01_miseq_nano/methylation_extraction/Histogram_BetaValues_PerSample_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
par(mfrow = c(2,4))
for(i in (seq(from = 1, to = 29, by = 4))){
hist(sub_all[,i], main = paste(gsub("_Methylation", "", colnames(sub_all)[i]), "CpG Methylation", sep = " : "), xlab = "Beta Value", ylab = "", axes = FALSE)
axis(1, cex = 0.8)
axis(2, las = 2, format(c(0,100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000), scientific = FALSE))
}
dev.off()

### calc average methylation per site
meth.mean<-apply(sub_all[,seq(from = 1, to = 29, by = 4)], 1, mean)
meth.median<-apply(sub_all[,seq(from = 1, to = 29, by = 4)], 1, median)
names(meth.mean)<-rownames(sub_all)
names(meth.median)<-rownames(sub_all)


locs<-unlist(strsplit(rownames(sub_all), "_"))
locs<-matrix(locs, ncol = 2, byrow = TRUE)
rownames(locs)<-rownames(sub_all)


pdf("/mnt/data1/isabel/RRBS/01_miseq_nano/methylation_extraction/Boxplot_MeanBetaValues_PerChr_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
boxplot(list(meth.mean[which(locs[,1] == "chr1")], meth.mean[which(locs[,1] == "chr2")], meth.mean[which(locs[,1] == "chr3")], meth.mean[which(locs[,1] == "chr4")],
meth.mean[which(locs[,1] == "chr5")],meth.mean[which(locs[,1] == "chr6")],meth.mean[which(locs[,1] == "chr7")],meth.mean[which(locs[,1] == "chr8")],
meth.mean[which(locs[,1] == "chr9")],meth.mean[which(locs[,1] == "chr10")],meth.mean[which(locs[,1] == "chr11")],meth.mean[which(locs[,1] == "chr12")], 
meth.mean[which(locs[,1] == "chr13")],meth.mean[which(locs[,1] == "chr14")],meth.mean[which(locs[,1] == "chr15")],meth.mean[which(locs[,1] == "chr16")],
meth.mean[which(locs[,1] == "chr17")],meth.mean[which(locs[,1] == "chr18")],meth.mean[which(locs[,1] == "chr19")],meth.mean[which(locs[,1] == "chr20")],
meth.mean[which(locs[,1] == "chrX")],meth.mean[which(locs[,1] == "chrM")]), names = c(1:20, "X", "M"), ylab = "Mean Beta Value", main = "CpG Methylation")

dev.off()

pdf("/mnt/data1/isabel/RRBS/01_miseq_nano/methylation_extraction/Boxplot_MedianBetaValues_PerChr_CpGMethylation_minReadDepth10.pdf", paper = "a4r", pointsize = 10, width = 0, height = 0)
boxplot(list(meth.median[which(locs[,1] == "chr1")], meth.median[which(locs[,1] == "chr2")], meth.median[which(locs[,1] == "chr3")], meth.median[which(locs[,1] == "chr4")],
meth.median[which(locs[,1] == "chr5")],meth.median[which(locs[,1] == "chr6")],meth.median[which(locs[,1] == "chr7")],meth.median[which(locs[,1] == "chr8")],
meth.median[which(locs[,1] == "chr9")],meth.median[which(locs[,1] == "chr10")],meth.median[which(locs[,1] == "chr11")],meth.median[which(locs[,1] == "chr12")], 
meth.median[which(locs[,1] == "chr13")],meth.median[which(locs[,1] == "chr14")],meth.median[which(locs[,1] == "chr15")],meth.median[which(locs[,1] == "chr16")],
meth.median[which(locs[,1] == "chr17")],meth.median[which(locs[,1] == "chr18")],meth.median[which(locs[,1] == "chr19")],meth.median[which(locs[,1] == "chr20")],
meth.median[which(locs[,1] == "chrX")],meth.median[which(locs[,1] == "chrM")]), names = c(1:20, "X", "M"), ylab = "Median Beta Value", main = "NonCpG Methylation")

dev.off()


all.herit<-read.csv("MethylationData_CpG_PresentAllSamples.csv", stringsAsFactors = FALSE, row.names = 1)

### denisty plots of output
density_meancovA17 <- density(all.herit[,2], bw = 0.05, from = 0, to = 1)
density_meancovA18 <- density(all.herit[,6], bw = 0.05, from = 0, to = 1)
density_meancovA19 <- density(all.herit[,10], bw = 0.05, from = 0, to = 1)
density_meancovA20 <- density(all.herit[,14], bw = 0.05, from = 0, to = 1)
density_meancovA21 <- density(all.herit[,18], bw = 0.05, from = 0, to = 1)
density_meancovA22 <- density(all.herit[,22], bw = 0.05, from = 0, to = 1)
density_meancovA23 <- density(all.herit[,26], bw = 0.05, from = 0, to = 1)
density_meancovA24 <- density(all.herit[,30], bw = 0.05, from = 0, to = 1)

jpeg('coverageA17.jpg')
plot(density_meancovA17, xlab = "Methylation A17", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA18.jpg')
plot(density_meancovA18, xlab = "Methylation A18", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA19.jpg')
plot(density_meancovA19, xlab = "Methylation A19", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA20.jpg')
plot(density_meancovA20, xlab = "Methylation A20", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA21.jpg')
plot(density_meancovA21, xlab = "Methylation A21", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA22.jpg')
plot(density_meancovA22, xlab = "Methylation A22", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA23.jpg')
plot(density_meancovA23, xlab = "Methylation A23", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()

jpeg('coverageA24.jpg')
plot(density_meancovA24, xlab = "Methylation A24", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")
dev.off()



####################### OTHER NOTES / TRIES #######################



source("https://bioconductor.org/biocLite.R")
biocLite("bsseq")
library(bsseq)

sample_names <- c("A17", "A24", "A26", "A28", "B29", "B30", "B31", "B32")

infile <- list.files(".", pattern = "^.*_r1_trimmed_bismark_bt2.bismark.cov.gz", all.files=TRUE, full.names=TRUE, recursive=FALSE, ignore.case=FALSE, include.dirs=FALSE, no..=TRUE)
bismarkBSseq <- read.bismark(files=infile, sampleNames=sample_names, rmZeroCov=FALSE, strandCollapse=FALSE, fileType=c("cov", "oldBedGraph", "cytosineReport"), mc.cores=1, verbose=TRUE)
bismarkBSseq # An object of type 'BSseq'
# ? a (matrix) of Cov (Coverage) values, describing the total number of reads covering a single loci. Each row in this matrix is a methylation loci and each column is a sample.

colnames(bismarkBSseq)

## https://www.bioconductor.org/packages/3.7/bioc/vignettes/bsseq/inst/doc/bsseq.html

# The genomic positions are stored as a GRanges object GRanges are general genomic regions
# Objects of class BSseq contain a GRanges object with the genomic locations.
# This GRanges object can be obtained by granges()
# For example, the first 4 loci in the Lister data are:
head(granges(bismarkBSseq), n = 4)

start(bismarkBSseq)
end(bismarkBSseq)
seqnames(bismarkBSseq) # chromosome names

# These objects also contains a phenoData object for sample pheno data.
sampleNames(bismarkBSseq)
pData(bismarkBSseq)

dim(bismarkBSseq)
ncol(bismarkBSseq) # number of columns; number of samples
nrow(bismarkBSseq) # number of rows; number of methylation loci

# BSseq() instantiates an object of class BSseq
# Genomic locations are passed in, either as a  GRanges object (argument gr) or as chromosome and location vectors (arguments chr and  pos).
# The arguments M and Cov accepts matrices, and it is possible to directly give it a  phenoData object.

# Obtaining coverage (methylation)
# Coverage, either Cov or M values, are obtained by getCoverage(), using the type argument:
head(getCoverage(bismarkBSseq, type = "Cov"))
head(getCoverage(bismarkBSseq, type = "M"))

# Methylation estimates can be obtained in the same way, using the function getMeth()
# If type is set to raw the function returns simple single-loci methylation estimates (which are M/Cov).
# To get smoothed estimates, the BSseq object needs to have been smoothed using Bsmooth, and  type set to smooth (default).

bismarkBSseq_ordered <- orderBSseq(bismarkBSseq)
granges(bismarkBSseq_ordered)

bismarkBSseq_smooth <- BSmooth(bismarkBSseq_ordered)
head(granges(bismarkBSseq_smooth))
start(bismarkBSseq_smooth)
end(bismarkBSseq_smooth)
seqnames(bismarkBSseq_smooth)
sampleNames(bismarkBSseq_smooth)
pData(bismarkBSseq_smooth)
dim(bismarkBSseq_smooth)
ncol(bismarkBSseq_smooth)
nrow(bismarkBSseq_smooth)
head(getCoverage(bismarkBSseq_smooth, type = "Cov"))
head(getCoverage(bismarkBSseq_smooth, type = "M"))

plotRegion(bismarkBSseq_smooth, col = c("black", "blue", "green", "red", "yellow", "pink", "orange", "purple")) # NO IDEA HOW TO CHANGE WIDTH!!!!

plotManyRegions(bismarkBSseq_smooth, col = c("black", "blue", "green", "red", "yellow", "pink", "orange", "purple"))


# RnBeads
# http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf
data.dir <- "/mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/" # Directory where your data is located




# read.bsmooth()
dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE) # "." means current directory
read.bsmooth(dirs, sampleNames = NULL, seqnames = NULL, returnRaw = FALSE, qualityCutoff = 20, rmZeroCov = FALSE, verbose = TRUE)


# An analysis typically consists of the following steps:
# 1. Smoothing, using the function BSmooth().
# 2. Compute t-statistics, using the function BSmooth.tstat(). This converts the BSseq object into a BSseqTstat object.
# 3. Threshold these t-statistics to identify DMRs, using the function dmrFinder(), returning a simple data.frame.

# It is usually a good idea to look at the smoothed data either before or after identifying DMRs. This can be done using the functions plotRegion() and plotManyRegions()

















## RnBeads (http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf)
# RnBeads also supports mouse and rat methylome analysis through its companion packages (RnBeads.mm9,RnBeads.mm10,RnBeads.rn5)
# Loading is performed with bed file input
# RnBeads expects sequencing data to be bed-formatted or have a similar coordinate file format
biocLite("RnBeads.mm10")
library(RnBeads.mm10)


### Running the modules step by step
# Directory where your data is located
data.dir <- "/mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/bed/"
bed.dir <- file.path(data.dir, "bs.bed.dir")
sample.annotation <- file.path(data.dir, "samples.csv")

# Directory where the output should be written to
analysis.dir <- "/mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/bed/"
# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports_details")
rnb.initialize.reports(report.dir)
# Set some analysis options
#rnb.options(filtering.sex.chromosomes.removal = TRUE, identifiers.column = "Sample_ID")
## Restrict logging to the console only
#logger.start(fname = NA)

## Data import
data.source <- "/mnt/data1/isabel/RRBS/01_miseq_nano/trim_galore/bed/"
result <- rnb.run.import(data.source=data.source, data.type="data.dir", dir.reports=report.dir)
rnb.set <- result$rnb.set
## Quality Control
rnb.run.qc(rnb.set, report.dir)

## Preprocessing
rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports=report.dir)$rnb.set

## Data export
rnb.options(export.to.csv = TRUE)
rnb.run.tnt(rnb.set, report.dir)

## Exploratory analysis
rnb.run.exploratory(rnb.set, report.dir)

## Differential methylation
rnb.run.differential(rnb.set, report.dir)




data.dir <- "~/RnBeads/data/Ziller2011_PLoSGen_450K"
     idat.dir <- file.path(data.dir, "idat")
     sample.annotation <- file.path(data.dir, "sample_annotation.csv")
     data.source <- c(idat.dir, sample.annotation)
     rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")





BS_A17 <- bismarkBSseq_smooth[,1]
save(BS_A17, file = "BS_A17.csv")
BS_A18 <- bismarkBSseq_smooth[,2]
save(BS_A18, file = "BS_A18.csv")
BS_A19 <- bismarkBSseq_smooth[,3]
save(BS_A19, file = "BS_A19.csv")
BS_A20 <- bismarkBSseq_smooth[,4]
save(BS_A20, file = "BS_A20.csv")
BS_A21 <- bismarkBSseq_smooth[,5]
save(BS_A21, file = "BS_A21.csv")
BS_A22 <- bismarkBSseq_smooth[,6]
save(BS_A22, file = "BS_A22.csv")
BS_A23 <- bismarkBSseq_smooth[,7]
save(BS_A23, file = "BS_A23.csv")
BS_A24 <- bismarkBSseq_smooth[,8]
save(BS_A24, file = "BS_A24.csv")

A17 <- read.csv("Heritability/HeritabilityEstimatesFromACMModel_All.csv", stringsAsFactors = FALSE, row.names = 1)

density_A17 <- density(BS_A17[,1], bw = 0.05, from = 0, to = 1)
par(mfrow = c(1,3))
par(mar = c(5,5,0.5,1))
plot(density_A17, xlab = "A", cex.axis = 1.5, cex.lab = 1.5, main = "", xaxs = "i")



## Using objects of class BSseq


getCoverage(BSseq=bismarkBSseq, regions=NULL, type=c("Cov", "M"), what=c("perBase", "perRegionAverage", "perRegionTotal"))





