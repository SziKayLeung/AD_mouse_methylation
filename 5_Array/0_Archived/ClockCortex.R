### Similiar to Clock.R script however only the DNAmAgeClockCortex will be used in this script


#### Does Aging in TG differ to WT

library(tidyr)
library(dplyr)
library(lme4)
library(mltools)
#################### ECX rTg4510
## load data

## read data in
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame.rdat") #load this in to have passed QC samples only as samplesheet has failed samples
pheno <- QCmetrics[which(QCmetrics$Tissue %in% c("CortexEntorhinalis", "Hippocampus")),] #take tissue samples further


samplesheet <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Horvath/MammalianArrayNormalizationTools/Clocks/OutputMouseClocks.csv", header = T,stringsAsFactors = F)
samplesheet <- samplesheet[,c(1:50)] #filter unnecessary values
samplesheet <- samplesheet[which(samplesheet$Basename %in% pheno$Basename),] #filter to just tissue samples to join with pathology data
samplesheet$SampleID <- gsub(".*_", "",samplesheet$ExternalSampleID)
## Add pathology data in
J20_pathdata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
rTg4510_pathdata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)

#Pathology knows which samples are paired so add the column for it now
for(i in 1:nrow(J20_pathdata)){
  J20_pathdata[i, "MouseID"] <- paste("J20_M", i, sep="")
}

for(i in 1:nrow(rTg4510_pathdata)){
  rTg4510_pathdata[i, "MouseID"] <- paste("rTg4510_M", i, sep="")
}

#split the tissue and model into sep data
J20_ecx <- J20_pathdata[,c("Sample_ID_ECX","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_ECX", "MouseID" )]
colnames(J20_ecx) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")
J20_hip <- J20_pathdata[,c("Sample_ID_HIP","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_HIP", "MouseID")]
colnames(J20_hip) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")

rTg4510_ecx <- rTg4510_pathdata[,c("Sample_ID_ECX","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_ECX", "MouseID")]
colnames(rTg4510_ecx) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")
rTg4510_hip <- rTg4510_pathdata[,c("Sample_ID_HIP","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_HIP", "MouseID")]
colnames(rTg4510_hip) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")

pathology <- rbind(J20_ecx, J20_hip, rTg4510_ecx, rTg4510_hip)
samplesheet <- merge(samplesheet, pathology, by = "SampleID")
samplesheet <- samplesheet[,-53]
colnames(samplesheet)[colnames(samplesheet) == 'Genotype.x'] <- 'Genotype'

samplesheet <- separate(data = samplesheet, col = Basename, 
                        into = c("Chip", "Chip_Position"), sep="_")
samplesheet$Genotype <- factor(samplesheet$Genotype, levels = c("WT","TG"))

#create a results table
res <- matrix(NA, nrow = 12, ncol= 4)
colnames(res) <- c("rTg4510_ECX", "rTg4510_HIP", "J20_ECX", "J20_HIP")
rownames(res) <- c("Genotype.Intercept", "Genotype.Beta", "Genotype.SE", "Genotype.P",
                   "Genotype.Age.Intercept","Genotype.Age.Beta","Genotype.Age.SE", "Genotype.Age.P", 
                   "Pathology.Intercept","Pathology.Beta", "Pathology.SE", "Pathology.P")


############### Functions for plots as it is repetitive and keep mislabelling 


#Plots for WT samples
plotWTSamples <- function(dat, clocks, tissue, model){
  dat2 <- dat[which(dat$Genotype == "WT"),]
  par(mfrow = c(2,3))
  for(clock in clocks){
    corr <- cor(dat2$Age, dat2[, clock])
    err <- rmse(preds = dat2[, clock], actuals = dat2$Age)
    plot(dat2$Age, dat2[, clock], col = dat2$GenotypeCol,
         xlab = "Age", ylab = as.character(clock),
         main = paste(paste(model, tissue, sep = " "), clock, sep = " - "))
    mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
    mtext(paste("rmse = ", signif(err,3)), side = 3, adj = 0)
    abline(0,1)
  }
}

#Plots for TGs and WTs samples
plotSamples <- function(dat, clocks, tissue, model){
  par(mfrow = c(2,3))
  for(clock in clocks){
    corr <- cor(dat$Age, dat[, clock])
    err <- rmse(preds = dat[, clock], actuals = dat$Age)
    plot(dat$Age, dat[, clock], col = dat$GenotypeCol,
         xlab = "Age", ylab = as.character(clock),
         main = paste(paste(model, tissue, sep = " "), clock, sep = " - "))
    mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
    mtext(paste("rmse = ", signif(err,3)), side = 3, adj = 0)
    abline(0,1)
  }
}



pdf("DNAmAgeClockCortex.pdf", w = 10)
########################################### filter for ecx rTg4510
Tissue <- "CortexEntorhinalis" 
Model  <- "rTg4510"
TG_colour <- "#00AEC9"

dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]
dat$Genotype <- factor(dat$Genotype, levels = c("WT","TG"))
dat$GenotypeCol <-  ifelse(dat$Genotype == "WT", "black", TG_colour)

## make a plot of actual age and DNAmAge for each type
clocks <- c("DNAmAgeClockCortex", "DNAmAgeClockBrain" , 
            "DNAmAgeClockPan","DNAmAgeClockLiver",
            "DNAmAgeClockBlood","DNAmAgeClockHeart")

## WT samples only
plotWTSamples(dat, clocks, tissue = "ECX", model = "rTg4510")

## all samples WT and TGs
plotSamples(dat, clocks, tissue = "ECX", model = "rTg4510")



## DNAmAge and terms
par(mfrow = c(1,3))
plot(dat$Genotype, dat$DNAmAgeClockCortex,
     xlab = "Genotype", ylab = "DNAmAgeClockCortex")
points(dat$Genotype, dat$DNAmAgeClockCortex, col = dat$GenotypeCol)

ages <- sort(unique(dat$Age_months))
means_WT <- c()
means_TG <- c()
for (age in unique(dat$Age_months)) {
  means_WT <- c(means_WT, mean(dat[which(dat$Age_months==age & dat$Genotype=="WT"),"DNAmAgeClockCortex"]))
  means_TG <- c(means_TG, mean(dat[which(dat$Age_months==age & dat$Genotype=="TG"),"DNAmAgeClockCortex"]))
}
plot(dat$Age_months, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Age months", ylab = "DNAmAgeClockCortex")
lines(ages, means_WT, lty=2, col="black")
lines(ages, means_TG, lty=2, col=TG_colour)

plot(dat$Pathology, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Pathology", ylab = "DNAmAgeClockCortex")
abline(lm(DNAmAgeClockCortex ~  Pathology, data=dat))

## Now run the linear models
lmmodel <- lm(DNAmAgeClockCortex ~  Genotype, data=dat)
summary(lmmodel)$coef
res[1,1] <- summary(lmmodel)$coef["(Intercept)",1]
res[2,1] <- summary(lmmodel)$coef["GenotypeTG",1]
res[3,1] <- summary(lmmodel)$coef["GenotypeTG",2]
res[4,1] <- summary(lmmodel)$coef["GenotypeTG",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Genotype*Age, data=dat)
summary(lmmodel)$coef
res[5,1] <- summary(lmmodel)$coef["(Intercept)",1]
res[6,1] <- summary(lmmodel)$coef["GenotypeTG:Age",1]
res[7,1] <- summary(lmmodel)$coef["GenotypeTG:Age",2]
res[8,1] <- summary(lmmodel)$coef["GenotypeTG:Age",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Pathology, data=dat)
summary(lmmodel)$coef
res[9,1] <- summary(lmmodel)$coef["(Intercept)",1]
res[10,1] <- summary(lmmodel)$coef["Pathology",1]
res[11,1] <- summary(lmmodel)$coef["Pathology",2]
res[12,1] <- summary(lmmodel)$coef["Pathology",4]

############################################ filter for hip rTg4510
Tissue <- "Hippocampus" 
Model  <- "rTg4510"
TG_colour <- "#00AEC9"

dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]
dat$Genotype <- factor(dat$Genotype, levels = c("WT","TG"))
dat$GenotypeCol <-  ifelse(dat$Genotype == "WT", "black", TG_colour)

## make a plot of actual age and DNAmAge for each type
clocks <- c("DNAmAgeClockCortex", "DNAmAgeClockBrain" , 
            "DNAmAgeClockPan","DNAmAgeClockLiver",
            "DNAmAgeClockBlood","DNAmAgeClockHeart")

## WT samples only
plotWTSamples(dat, clocks, tissue = "HIP", model = "rTg4510")

## all samples WT and TGs
plotSamples(dat, clocks, tissue = "HIP", model = "rTg4510")


## DNAmAge and terms
par(mfrow = c(1,3))
plot(dat$Genotype, dat$DNAmAgeClockCortex,
     xlab = "Genotype", ylab = "DNAmAgeClockCortex")
points(dat$Genotype, dat$DNAmAgeClockCortex, col = dat$GenotypeCol)

ages <- sort(unique(dat$Age_months))
means_WT <- c()
means_TG <- c()
for (age in unique(dat$Age_months)) {
  means_WT <- c(means_WT, mean(dat[which(dat$Age_months==age & dat$Genotype=="WT"),"DNAmAgeClockCortex"]))
  means_TG <- c(means_TG, mean(dat[which(dat$Age_months==age & dat$Genotype=="TG"),"DNAmAgeClockCortex"]))
}
plot(dat$Age_months, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Age months", ylab = "DNAmAgeClockCortex")
lines(ages, means_WT, lty=2, col="black")
lines(ages, means_TG, lty=2, col=TG_colour)

plot(dat$Pathology, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Pathology", ylab = "DNAmAgeClockCortex")
abline(lm(DNAmAgeClockCortex ~  Pathology, data=dat))

## Now run the linear models
lmmodel <- lm(DNAmAgeClockCortex ~  Genotype , data=dat)
summary(lmmodel)$coef
res[1,2] <- summary(lmmodel)$coef["(Intercept)",1]
res[2,2] <- summary(lmmodel)$coef["GenotypeTG",1]
res[3,2] <- summary(lmmodel)$coef["GenotypeTG",2]
res[4,2] <- summary(lmmodel)$coef["GenotypeTG",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Genotype*Age, data=dat)
summary(lmmodel)$coef
res[5,2] <- summary(lmmodel)$coef["(Intercept)",1]
res[6,2] <- summary(lmmodel)$coef["GenotypeTG:Age",1]
res[7,2] <- summary(lmmodel)$coef["GenotypeTG:Age",2]
res[8,2] <- summary(lmmodel)$coef["GenotypeTG:Age",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Pathology, data=dat)
summary(lmmodel)$coef
res[9,2] <- summary(lmmodel)$coef["(Intercept)",1]
res[10,2] <- summary(lmmodel)$coef["Pathology",1]
res[11,2] <- summary(lmmodel)$coef["Pathology",2]
res[12,2] <- summary(lmmodel)$coef["Pathology",4]

############################################ filter for ecx J20
Tissue <- "CortexEntorhinalis" 
Model  <- "J20"
TG_colour <- "#FF5A62"

dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]
dat$Genotype <- factor(dat$Genotype, levels = c("WT","TG"))
dat$GenotypeCol <-  ifelse(dat$Genotype == "WT", "black", TG_colour)

## make a plot of actual age and DNAmAge for each type
clocks <- c("DNAmAgeClockCortex", "DNAmAgeClockBrain" , 
            "DNAmAgeClockPan","DNAmAgeClockLiver",
            "DNAmAgeClockBlood","DNAmAgeClockHeart")
## WT samples only
plotWTSamples(dat, clocks, tissue = "ECX", model = "J20")

## all samples WT and TGs
plotSamples(dat, clocks, tissue = "ECX", model = "J20")

## DNAmAge and terms
par(mfrow = c(1,3))
plot(dat$Genotype, dat$DNAmAgeClockCortex,
     xlab = "Genotype", ylab = "DNAmAgeClockCortex")
points(dat$Genotype, dat$DNAmAgeClockCortex, col = dat$GenotypeCol)

ages <- sort(unique(dat$Age_months))
means_WT <- c()
means_TG <- c()
for (age in unique(dat$Age_months)) {
  means_WT <- c(means_WT, mean(dat[which(dat$Age_months==age & dat$Genotype=="WT"),"DNAmAgeClockCortex"]))
  means_TG <- c(means_TG, mean(dat[which(dat$Age_months==age & dat$Genotype=="TG"),"DNAmAgeClockCortex"]))
}
plot(dat$Age_months, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Age months", ylab = "DNAmAgeClockCortex")
lines(ages, means_WT, lty=2, col="black")
lines(ages, means_TG, lty=2, col=TG_colour)

plot(dat$Pathology, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Pathology", ylab = "DNAmAgeClockCortex")
abline(lm(DNAmAgeClockCortex ~  Pathology, data=dat))


## Now run the linear models
lmmodel <- lm(DNAmAgeClockCortex ~  Genotype, data=dat)
summary(lmmodel)$coef
res[1,3] <- summary(lmmodel)$coef["(Intercept)",1]
res[2,3] <- summary(lmmodel)$coef["GenotypeTG",1]
res[3,3] <- summary(lmmodel)$coef["GenotypeTG",2]
res[4,3] <- summary(lmmodel)$coef["GenotypeTG",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Genotype*Age, data=dat)
summary(lmmodel)$coef
res[5,3] <- summary(lmmodel)$coef["(Intercept)",1]
res[6,3] <- summary(lmmodel)$coef["GenotypeTG:Age",1]
res[7,3] <- summary(lmmodel)$coef["GenotypeTG:Age",2]
res[8,3] <- summary(lmmodel)$coef["GenotypeTG:Age",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Pathology, data=dat)
summary(lmmodel)$coef
res[9,3] <- summary(lmmodel)$coef["(Intercept)",1]
res[10,3] <- summary(lmmodel)$coef["Pathology",1]
res[11,3] <- summary(lmmodel)$coef["Pathology",2]
res[12,3] <- summary(lmmodel)$coef["Pathology",4]



############################################ filter for hip J20
Tissue <- "Hippocampus" 
Model  <- "J20"
TG_colour <- "#FF5A62"

dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]
dat$Genotype <- factor(dat$Genotype, levels = c("WT","TG"))
dat$GenotypeCol <-  ifelse(dat$Genotype == "WT", "black", TG_colour)

## make a plot of actual age and DNAmAge for each type
clocks <- c("DNAmAgeClockCortex", "DNAmAgeClockBrain" , 
            "DNAmAgeClockPan","DNAmAgeClockLiver",
            "DNAmAgeClockBlood","DNAmAgeClockHeart")

## WT samples only
plotWTSamples(dat, clocks, tissue = "HIP", model = "J20")

## all samples WT and TGs
plotSamples(dat, clocks, tissue = "HIP", model = "J20")

## DNAmAge and terms
par(mfrow = c(1,3))
plot(dat$Genotype, dat$DNAmAgeClockCortex,
     xlab = "Genotype", ylab = "DNAmAgeClockCortex")
points(dat$Genotype, dat$DNAmAgeClockCortex, col = dat$GenotypeCol)

ages <- sort(unique(dat$Age_months))
means_WT <- c()
means_TG <- c()
for (age in unique(dat$Age_months)) {
  means_WT <- c(means_WT, mean(dat[which(dat$Age_months==age & dat$Genotype=="WT"),"DNAmAgeClockCortex"]))
  means_TG <- c(means_TG, mean(dat[which(dat$Age_months==age & dat$Genotype=="TG"),"DNAmAgeClockCortex"]))
}
plot(dat$Age_months, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Age months", ylab = "DNAmAgeClockCortex")
lines(ages, means_WT, lty=2, col="black")
lines(ages, means_TG, lty=2, col=TG_colour)

plot(dat$Pathology, dat$DNAmAgeClockCortex, col = dat$GenotypeCol, 
     xlab = "Pathology", ylab = "DNAmAgeClockCortex")
abline(lm(DNAmAgeClockCortex ~  Pathology, data=dat))


## Now run the linear models
lmmodel <- lm(DNAmAgeClockCortex ~  Genotype, data=dat)
summary(lmmodel)$coef
res[1,4] <- summary(lmmodel)$coef["(Intercept)",1]
res[2,4] <- summary(lmmodel)$coef["GenotypeTG",1]
res[3,4] <- summary(lmmodel)$coef["GenotypeTG",2]
res[4,4] <- summary(lmmodel)$coef["GenotypeTG",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Genotype*Age, data=dat)
summary(lmmodel)$coef
res[5,4] <- summary(lmmodel)$coef["(Intercept)",1]
res[6,4] <- summary(lmmodel)$coef["GenotypeTG:Age",1]
res[7,4] <- summary(lmmodel)$coef["GenotypeTG:Age",2]
res[8,4] <- summary(lmmodel)$coef["GenotypeTG:Age",4]

lmmodel <- lm(DNAmAgeClockCortex ~  Pathology, data=dat)
summary(lmmodel)$coef
res[9,4] <- summary(lmmodel)$coef["(Intercept)",1]
res[10,4] <- summary(lmmodel)$coef["Pathology",1]
res[11,4] <- summary(lmmodel)$coef["Pathology",2]
res[12,4] <- summary(lmmodel)$coef["Pathology",4]

dev.off()


write.csv(res, "DNAmAgeCloclCortex_linearmodelresults.csv")

