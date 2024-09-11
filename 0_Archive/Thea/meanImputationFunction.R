## Mean imputation of Isabel's samples to make the mouse clock work on them

# Turn this into a function to be run within a shell script

## INPUT: a single .cov file 

## FUNCTION: Carry out mean imputation of the sample and keep only probes 
##           needed for the MouseEpigeneticClock

## OUTPUT: R matrix containing clock probes

## TO DO AFTER IN SHELL SCRIPT: save each matrix as .cov and run through clock

## Dorothea Seiler Vellame 21-08-2018


# load data for testing #############################################################
# setwd("/mnt/data1/Thea/EpigeneticClock/MouseEpigeneticClock/bismarkFiles")
# 
# # load single file to practice on
# temp=list.files(pattern="*.cov")
# sampleList<-list()
# for (i in temp[1]){
#   sampleList[[i]]<-assign(substr(i,6,8),
#                           as.data.frame(read.table(i,header=FALSE,sep="\t",
#                           stringsAsFactors = FALSE,quote="")))
# }
# covFile<- K17
# rm(K17,i,sampleList,temp)
#####################################################################################

meanImputationFunction <- function(covFile){
  
  # load the sites needed for prediction
  setwd("/mnt/data1/Thea/EpigeneticClock/MouseEpigeneticClock")
  RdataF = "./PredictionPackage_20170215.Rdata"
  load(RdataF)
  rm(betas,qnTarget,RdataF,rowMean,rowStDev)
  
  # find all prediction sites included in the sample
  includedInSites<-paste(substring(covFile[,1],4), covFile[,2],sep=":") %in% sitesForPrediction
  
  # take subset of covFile of included prediction sites
  covFileSub<-covFile[includedInSites,]
  
  # create list of missing prediction sites
  missingFromSample<-!(sitesForPrediction %in% paste(substring(covFileSub[,1],4), covFileSub[,2],sep=":"))
  
  # create bottom of covFileSub for missing sites
  covFilePlus<-matrix(ncol=6,nrow=sum(missingFromSample))
  for (i in 1:sum(missingFromSample)){
    covFilePlus[i,1]<-paste("chr",strsplit(sitesForPrediction[missingFromSample],":")[[i]][1],sep="")
    covFilePlus[i,2]<-strsplit(sitesForPrediction[missingFromSample],":")[[i]][2]
  }
  
  # find mean methylation value and fill covFilePlus
  methylationMean<-mean(covFile[,4])
  covFilePlus[,4]<-as.numeric(methylationMean)
  covFilePlus[,5:6]<-as.numeric(6)
  
  # bind the covFile to include all prediction sites
  covFileAll<-rbind(covFileSub,covFilePlus)
  
  # check that coverage is greater than 5 for each site
  for (i in 1:nrow(covFileAll)){
    if ((as.numeric(covFileAll[i,5])+as.numeric(covFileAll[i,6]))<5){
      covFileAll[i,4]<-as.numeric(methylationMean)
      covFileAll[i,5]<-as.numeric(6)
    }
  }
  
  class(covFileAll[,4])<-"numeric"
  class(covFileAll[,5])<-"numeric"
  class(covFileAll[,6])<-"numeric"
  
  return(covFileAll)
}









# setwd("/mnt/data1/Thea/EpigeneticClock/MouseEpigeneticClock/bismarkFiles")
# save(covFileAll,file="covFile.Rdata")
# 
# 
# # check that the file is okay
# rm(list=ls(all=TRUE))
# load(file)
# load("covFile.Rdata")
# head(covFileAll)
