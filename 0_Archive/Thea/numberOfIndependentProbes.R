## calculate number of unique probes for RRBS data
# do so post filtering so that you ignore probes with too low read depth

# Dorothea Seiler Vellame 26-06-2019

# load data
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("J20MatrixFiltered.Rdata")
J20 = RRBS
load("rTg4510MatrixFiltered.Rdata")
rTg = RRBS

# create list that has each chromosome as an object
# methylation cols only

chrom = unlist(strsplit(as.character(J20[,1]), split = "_"))
chrom = chrom[seq(1,length(chrom),2)]
J20chrom = list()
for (i in 1:length(levels(as.factor(chrom)))){
  J20chrom[[i]] = J20[which(chrom==levels(as.factor(chrom))[i]),2:64]
  rownames(J20chrom[[i]]) = J20[which(chrom==levels(as.factor(chrom))[i]),1]
}


chrom = unlist(strsplit(as.character(rTg[,1]), split = "_"))
chrom = chrom[seq(1,length(chrom),2)]
rTgchrom = list()
for (i in 1:length(levels(as.factor(chrom)))){
  rTgchrom[[i]] = rTg[which(chrom==levels(as.factor(chrom))[i]),2:63]
  rownames(rTgchrom[[i]]) = rTg[which(chrom==levels(as.factor(chrom))[i]),1]
}




## formulas from Moskvina paper
Moskvina<-function(chromData, alpha){
  # create correlation matrix from the data
  data = cor(t(chromData), method = "pearson", use = "p")
  
  # Moskvina method
  kj<-rep(NA, (nrow(data)-1))
  for(j in 1:(nrow(data)-1)){
    kj[j]<-sqrt(1-max(abs(data[(j+1),1:j]),na.rm = TRUE)^(-1.31*log10(alpha)))
  }
  return((1+sum(kj)))
}



uniqueJ20 = matrix(nrow = length(J20chrom), ncol = 2)
colnames(uniqueJ20) = c("unique CpGs", "Total CpGs")

for (i in 1:length(J20chrom)){
  uniqueJ20[i,2] = dim(J20chrom[[i]])[1]
  uniqueJ20[i,1] = Moskvina(J20chrom[[i]], 0.05)
}


uniquerTg = matrix(nrow = length(rTgchrom), ncol = 2)
colnames(uniquerTg) = c("unique CpGs", "Total CpGs")

for (i in 1:length(rTgchrom)){
  uniquerTg[i,2] = dim(rTgchrom[[i]])[1]
  uniquerTg[i,1] = Moskvina(rTgchrom[[i]], 0.05)
}

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(uniqueJ20, uniquerTg, file ="NumberOfUniqueProbesPerChrom.Rdata")


# load data and plot 
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("NumberOfUniqueProbesPerChrom.Rdata")

uniqueJ20 = uniqueJ20[1:19,]
uniquerTg = uniquerTg[1:19,]
plot(uniqueJ20[,1],uniqueJ20[,2])
plot(uniquerTg[,1],uniquerTg[,2])

proportions = data.frame(c(uniqueJ20[,1]/uniqueJ20[,2]*100, uniquerTg[,1]/uniquerTg[,2]*100))
colnames(proportions) = "prop"

library(ggplot2)
ggplot(data = proportions, aes(x = "", y = prop)) + 
  geom_boxplot() +
  xlab("") + 
  ylab("Proportion of CpGs that are unique, %") +
  theme_grey(base_size = 18)
  
  
## make single heatmap of chromosome 1

data = cor(t(J20chrom[[8]]), method = "pearson", use = "p")

corrplot::corrplot(data[1:100,1:100])





