## Remove non-vaiable probes, NAs and filter by read depth

## INPUT: RRBSmatrix (output from RRBSMatrix.R)
##        RDThreshold - minimum coverage per sample site
##        NAmax - maximum number of NAs allowed per sample MUST BE LESS THAN TOTAL
##        densityPlots (default FALSE) - if TRUE will plot density plots of methylation 
##        before and after removal of non variable probes
##        SDThreshold - filtering of non-variabke probes by standard deviation
##        covCols = FALSE, when TRUE coverage columns will be kept

## OUTPUT: RRBSfiltered 

# Dorothea Seiler Vellame 06-07-2018

# 18-10-2018 INPUT added: groupPheno - Array of phenotype that group is being split by. Length and order
# must be equal to that of the samples and they must be factors

# Number of NAMax will be changed to sampleMin, the minimum number of samples per group. This must be larger
# or equal to the smallest group size

# 23-10-2018 INPUT added: pdfFileName This names the density plot pdf if densityPlots=TRUE. .pdf should be included in the name.

# suggested inputs
#setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
#source("RRBSMatrix.R")
#RRBSmatrix<-RRBSMatrix("/mnt/data1/Thea/IsabelSamples/data/covFolder")
#setwd("/mnt/data1/Thea/IsabelSamples/data")
#load("J20MatrixUnfiltered.Rdata")
#RDThreshold<-20
#NAMax<-1

RRBSFiltering <- function (RRBSmatrix, RDThreshold, groupPheno, sampleMin,
                           densityPlots=FALSE,
                           pdfFileName="densityFilteringPlots.pdf",
                           SDThreshold = 5,
                           covCols=FALSE){
  
  numberOfSamples<-(ncol(RRBSmatrix)-1)/2
  methylationCols<-2:(1+numberOfSamples)
  coverageCols<-(2+numberOfSamples):ncol(RRBSmatrix)
  
  # check that groupPheno is a factor
  if (class(groupPheno)!="factor"){
    stop("groupPheno must be class factor")
  }
  
  # checks that length of groupPheno matches the number of cols in RRBSmatrix
  if (length(groupPheno)!=numberOfSamples){
    stop("Number of phenotypes not equal to number of samples in matrix")
  }
  
  # ensure that sample min is greater than the smallest sample group
  if (sum(sampleMin>summary(groupPheno))!=0){
    stop("sampleMin larger than smallest group sample size. Decrease sampleMin")
  }
  
  
  
  ## filter by read depth #######################################
  RRBSM<-RRBSmatrix[,methylationCols]
  RRBSC<-RRBSmatrix[,coverageCols]
  k<-RRBSC<RDThreshold
  RRBSM[k]<-NA
  RRBSC[k]<-NA
  RRBSmatrix<-cbind(RRBSmatrix[,1],RRBSM,RRBSC)
  

  ## remove rows with too many NA's #############################
  ## remove rows where minimum number of samples aren't present per group
  
  # return row index of rows to be removed 
  NAFilterIndex <- function(RRBSrow,levels,sampleMin,groupPheno){
    groupKeep<-matrix(ncol=length(levels),nrow=1)
    for (i in 1:length(levels)){
      group<-which(groupPheno==levels[i])
      if (length(group)-sum(is.na(RRBSrow[group]))>=sampleMin){
        groupKeep[i]<-TRUE
      }else{ 
        groupKeep[i]<- FALSE
      }
    }
    if (sum(groupKeep)==length(levels)){
      return("keep")
    }else{
      return(0)
    }
  }
  
  
  levels<-levels(groupPheno)
  k<-apply(RRBSmatrix[,methylationCols],1,NAFilterIndex,levels,sampleMin,groupPheno)
  RRBSmatrix<-RRBSmatrix[which(k=="keep"),]
  
  
  ## remove non variable probes by standard deviation ###########
  sdMatrix<-apply(RRBSmatrix[,methylationCols],1,sd,na.rm=TRUE)
  k<-which(sdMatrix>SDThreshold)
  
  if (densityPlots==TRUE){
    pdf(file=pdfFileName)
    for (i in methylationCols){
      par(mfrow=c(2,1))
      plot(density(RRBSmatrix[,i],na.rm=TRUE), main=paste("methylation distribution in",colnames(RRBSmatrix)[i],"raw"))
      plot(density(unlist(RRBSmatrix[k,i]),na.rm=TRUE), main=paste("methylation distribution in",colnames(RRBSmatrix)[i],"filtered"))
    }
    dev.off()
  }
  RRBSmatrix<-RRBSmatrix[k,]
  
  # remove the coverage values if covCols=FALSE as they will no longer be needed
  if (covCols==TRUE){
    return(RRBSmatrix)
  }
  RRBSmatrix<-RRBSmatrix[,1:(numberOfSamples+1)]
  return(RRBSmatrix)
  
}