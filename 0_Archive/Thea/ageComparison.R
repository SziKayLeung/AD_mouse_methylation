data<-c("K17", 19.4,
"K18", 6.7,
"K19", -1.1,
"K20", 2.3,
"K21", 0.3,
"K23", 4.2,
"K24", 65.4,
"L17", 3.3,
"L18", 8.2,
"L19", 0,
"L21", 33.8,
"L22", 3.9,
"L23", -0.1,
"L24", 22.6,
"M17", 17.1,
"M18", 1.7,
"M19", 6.3,
"M20", 3.9,
"M21", 2.4,
"M22", -0.4,
"M23", 8.1,
"M24", 0.1,
"N24", 1.9,
"O18", 34.9,
"O19", 2.7,
"O20", 2.7,
"O21", 29.4,
"O22", 15.9,
"O23", -1.5,
"O24", 29.9,
"P17", -1,
"P18", -0.1,
"P19", 1.1,
"P20", 3.2,
"P21", 4.8,
"P22", 30,
"P23", 107.6,
"P24", 2.3,
"Q17", 19.6,
"Q18", 0,
"Q19", 2,
"Q20", 9.2,
"Q21", 2,
"Q22", 4,
"Q23", 1.2,
"Q24", 9,
"S17", 14.6,
"S18", 7.8,
"S19", 1,
"S20", 10.9,
"S21", 23.3,
"S22", 0.4,
"S23", 27.9,
"S24", 8.7,
"T17", 2.3,
"T18", 4.5,
"T19", -1.6,
"T20", -0.3,
"T21", 13.9,
"T22", 15.3,
"T23", -1.1,
"T24", 22.2)

names<-data[seq(1,123,2)]
predictedAge<-data[seq(2,124,2)]
pred<-cbind(names,predictedAge,matrix(nrow=62,ncol=3))

# load phenotypes
setwd("/mnt/data1/Thea/IsabelSamples/data")
phenotype<-read.csv("Tg4510_samples_table.csv")

pred[,3]<-phenotype[match(pred[,1],phenotype[,1]),5]
pred[,4]<-as.character(phenotype[match(pred[,1],phenotype[,1]),4])
pred[,5]<-as.numeric(pred[,3])*4 # convert months to weeks

pred<-pred[!is.na(pred[,3]),]

plot(pred[,5],pred[,2], col=ifelse(pred[,4]=="WT","blue","red"),
     xlab="Age (weeks)", pch=20, cex=2,cex.lab=1.3,cex.axis=1.5,
     ylab="Predicted age (weeks)")
legend("topright",unique(pred[,4]),col=c("blue","red"),pch=20,cex = 1.5)
cor(as.numeric(pred[,2]),as.numeric(pred[,3]))
# 0.01577138
# With mean imputation not a strong correlation! BUT you can't know how good the age predictor it self is!