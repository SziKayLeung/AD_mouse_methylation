setwd("/mnt/data1/Thea/IsabelSamples/data/")

load("J20AgeStats.Rdata")
load("J20GenotypeStats.Rdata")

J20Age<-RRBSAgeStat
J20Genotype<-RRBSGenotypeStat

load("rTg4510AgeStats.Rdata")
load("rTg4510GenotypeStats.Rdata")

rTg4510Age<-RRBSAgeStat
rTg4510Genotype<-RRBSGenotypeStat

## Check correlation between p-value and total n

plot(J20Age[,6],J20Age[,4])
plot(J20Genotype[,6],J20Genotype[,4])
plot(rTg4510Age[,6],rTg4510Age[,4])
plot(rTg4510Genotype[,6],rTg4510Genotype[,4])

## No relationship which is good!Actually, of the top 10000 have higher n which is the way you'd want it