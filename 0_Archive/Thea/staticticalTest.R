# Statistical analysis of RRBS data

## INPUT: RRBSm - matrix of RRBS data with sites as first column and each subsequent column as a sample
##        groupPheno - phenotype with col 1 genotype, col 2 age

## OUTPUT: outList - list holding 3 objects, containting the stats results for genotype, age and the
##                   interaction respectively


# Dorothea Seiler Vellame 06-11-2018  

# # load data and subset for testing
# setwd("/mnt/data1/Thea/IsabelSamples/data/")
# load("J20MatrixFiltered.Rdata")
# load("J20MatrixUnfiltered.Rdata")
# RRBSm<-RRBS[1:50,1:64]
# 
# #make smaller subset
# RRBSmBIG<-RRBSm
# RRBSm<-RRBSmBIG[,]
# 
# rm(RRBS)
# 
# phenotype<-read.csv("J20_phenotype_RRBS.csv")
# groupPheno<-phenotype[,c(4,5)]
# rm(phenotype)

statisticalTest<-function(RRBSm,groupPheno){
  
  # create list to hold output
  genotype<-matrix(nrow=nrow(RRBSm),ncol=5)
  age<-matrix(nrow=nrow(RRBSm),ncol=5)
  ageGenoInteraction<-matrix(nrow=nrow(RRBSm),ncol=5)
  
  outList<-list(genotype,age,ageGenoInteraction)
  
  
  # create function to be used in apply across all rows
  
  anovaExtract<-function(RRBSCol,groupPheno,RRBSnrow){
    
    dataAndPheno<-cbind(RRBSCol[2:length(RRBSCol)],groupPheno)
    
    fit<-lm(as.numeric(dataAndPheno[,1])~dataAndPheno[,2]+as.factor(dataAndPheno[,3])+dataAndPheno[,2]*as.factor(dataAndPheno[,3]))
    x<-anova(fit)
    
    returnMatrix<-matrix(ncol=5,nrow=3)
    returnMatrix[,1]<-x$`F value`[1:3]
    returnMatrix[,2]<-x$Df[1:3]
    returnMatrix[,4]<-x$`Pr(>F)`[1:3]
    returnMatrix[,5]<-p.adjust(returnMatrix[,4],method="fdr",n=RRBSnrow)                             
    returnMatrix[1,3]<-sqrt(x$`Mean Sq`[1]/x$`Mean Sq`[4])
    returnMatrix[2,3]<-sqrt(x$`Mean Sq`[2]/x$`Mean Sq`[4])
    returnMatrix[3,3]<-sqrt(x$`Mean Sq`[3]/x$`Mean Sq`[4])
    
    return(t(returnMatrix))
  }
  
  flatAnova<-apply(RRBSm,1,anovaExtract,groupPheno,nrow(RRBSm))
  rownames(flatAnova)<-c("F value","Degree of Freedom","Effect size","p value","FDR p","F value","Degree of Freedom","Effect size","p value","FDR p","F value","Degree of Freedom","Effect size","p value","FDR p")
  colnames(flatAnova)<-RRBSm[,1]
  
  
  
  # fill outlist from the apply output
  outList[[1]]<-t(flatAnova[1:5,])
  outList[[2]]<-t(flatAnova[6:10,])
  outList[[3]]<-t(flatAnova[11:15,])
  
  # sort each list by p-value (too many FDR's are 1)
  outList[[1]]<-outList[[1]][order(outList[[1]][,4]),]
  outList[[2]]<-outList[[2]][order(outList[[2]][,4]),]
  outList[[3]]<-outList[[3]][order(outList[[3]][,4]),]
  
  return(outList)
}
