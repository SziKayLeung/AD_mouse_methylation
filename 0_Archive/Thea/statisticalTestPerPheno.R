## A statistical function that uses only one genotype at a time

# Dorothea Seiler Vellame 07-01-2019

# testing data
# setwd("/mnt/data1/Thea/IsabelSamples/data")
# load("J20AgeMatrixFiltered.Rdata")
# 
# RRBSm<-RRBS[1:100,1:64]
# 
# # load the phenotype data and extract the groups column
# setwd("/mnt/data1/Thea/IsabelSamples/data/")
# phenotype<-read.csv("J20_phenotype_RRBS.csv")
# groupPheno<-as.factor(phenotype[,5]) 
# 
# RRBSRow<-RRBSm[1,]


statisticalTestPerPheno<-function(RRBSm,groupPheno){
  
  anovaExtract<- function(RRBSRow, groupPheno,RRBSnrow){
    
    dataAndPheno<-rbind(as.numeric(RRBSRow[2:length(RRBSRow)]),as.factor(groupPheno))  # create a data frame containing the wanted row with genotype data
    if (length(which(is.na(dataAndPheno[1,])))>0){                                # remove NA entries so that lm works
      dataAndPheno<-dataAndPheno[,-which(is.na(dataAndPheno[1,]))]
    }

    fit<-lm(as.numeric(dataAndPheno[1,])~as.factor(dataAndPheno[2,]))             # apply standard model
    
    x<-anova(fit)
    
    returnMatrix<-matrix(ncol=6,nrow=1)                                           # create return matrix and calculate wanted entries
    returnMatrix[1]<-x$`F value`[1]   
    returnMatrix[2]<-x$Df[1]
    returnMatrix[3]<-sqrt(x$`Mean Sq`[1]/x$`Mean Sq`[2])                          # absolute difference not included as it doesn't make sense for a phenotype with >2 groups
    returnMatrix[4]<-x$`Pr(>F)`[1]
    returnMatrix[5]<-p.adjust(returnMatrix[,4],method="fdr",n=RRBSnrow)                             
    returnMatrix[6]<-ncol(dataAndPheno)
    
    return(returnMatrix)
  }
  
  anovaOut<-as.data.frame(apply(RRBSm,1,anovaExtract,groupPheno,nrow(RRBSm)))            # configure output ot desired format
  rownames(anovaOut)<-c("F_value","Degree_of_Freedom","Effect_size","p_value","FDR_p", "n_Total")
  colnames(anovaOut)<-RRBSm[,1]
  anovaOut<-t(anovaOut)
  anovaOut[,c(2,6)]<-as.integer(anovaOut[,c(2,6)])
  anovaOut<-anovaOut[order(anovaOut[,4]),]                             # sort by p value
 
  return(anovaOut)
}
