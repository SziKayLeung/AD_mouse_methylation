# Isabel Castanho

# Raw reads: Total Sequences in fastqc.html
## J20
# 02_Hiseq_0144
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/02_Hiseq_0144/2605_*_r1_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*_r1.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/fastqc_02_Hiseq_0144.txt
# 03_HiSeq_0145
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/03_HiSeq_0145/2605_*_r1_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*_r1.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/fastqc_03_HiSeq_0145.txt

## rT4510
# 04_Hiseq_0146
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/04_Hiseq_0146/2605_*_r1_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*_r1.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/fastqc_04_Hiseq_0146.txt
# 05_Hiseq_0147
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/05_Hiseq_0147/2605_*_r1_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*_r1.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/fastqc_05_Hiseq_0147.txt


# Raw reads: Total Sequences in fastqc.html
## J20
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/J20/2605_*_r1_trimmed_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*_r1_trimmed.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/trimmed_fastqc_J20.txt
## rT4510
cat /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/Tg4510/2605_*_r1_trimmed_fastqc.html | egrep -oh '<td>Filename</td><td>2605_(.)*r1_trimmed.fq.gz</td>|<td>Total Sequences</td><td>[0-9]*</td>' > /mnt/data1/isabel/RRBS/RRBS-ADmice-ECX/trimmed_fastqc_Tg4510.txt


# t-test in R
R

setwd("/mnt/data1/isabel/RRBS/qc/")

# Tg4510
rawReads <- read.csv("/mnt/data1/isabel/RRBS/qc/Tg4510NumberRawReads.csv", row.names=1, stringsAsFactors=TRUE)
rawReads <- as.data.frame(rawReads)
colnames(rawReads) <- c("Genotype", "raw_reads")
t.test(rawReads$raw_reads ~ rawReads$Genotype) # where y is numeric and x is a binary factor # t = 0.11578, df = 54.262, p-value = 0.9083


# J20
rawReads <- read.csv("/mnt/data1/isabel/RRBS/qc/J20NumberRawReads.csv", row.names=1, stringsAsFactors=TRUE)
rawReads <- as.data.frame(rawReads)
colnames(rawReads) <- c("Genotype", "raw_reads")
t.test(rawReads$raw_reads ~ rawReads$Genotype) # where y is numeric and x is a binary factor # t = -0.76256, df = 59.418, p-value = 0.4487