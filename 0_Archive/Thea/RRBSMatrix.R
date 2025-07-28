## Make a table of methylation and coverage values from Isabel's data

## INPUT: file path to unzippped cov files (upload from harddrive)

## OUTPUT: RRBS matrix containing chromosome_startPosition, methylation 
##         value per sample and coverage value per sample. 
##         methylationCols and coverageCols containing the rownumbers of 
##         methylation and coverage data respectively.

# Dorothea Seiler Vellame 05-07-2018

# input file path:
# filePath="/mnt/data1/Thea/IsabelSamples/data/covFolder/"

RRBSMatrix<- function (filePath){
  setwd(filePath)
  temp=list.files(pattern="*.cov")
  
  ## read data into a list
  sampleList<-list()
  for (i in temp){
    sampleList[[i]]<-assign(substr(i,6,8),
                            as.data.frame(read.table(i,header=FALSE,sep="\t",
                                                     stringsAsFactors = FALSE,quote=""))) 
  }
  
  names<-substr(temp,6,8)
  numberOfSamples<-length(names)
  
  ## For each sample, format so that you have 3 cols: chr_startPosition, %m, cov
  # To use, loop through samples and names together so that they match
  RRBSColSelect<- function(RRBS, name){
    RRBStemp<-as.data.frame(matrix(ncol=3,nrow=nrow(RRBS)))
    RRBStemp[,1]<-paste(RRBS[,1],RRBS[,2],sep="_")
    RRBStemp[,2]<-RRBS[,4]
    RRBStemp[,3]<-RRBS[,5]+RRBS[,6]
    colnames(RRBStemp)<-c(paste("chr_start",name,sep="_"),paste("%m",name,sep="_"),
                          paste("cov",name,sep="_"))
    return(RRBStemp)
  }
  
  # loop through the sampleList applying ColSelect
  for(i in 1:numberOfSamples){
    sampleList[[i]]<-RRBSColSelect(sampleList[[i]],names[i])
  }
  
  ## create unique list of all sites
  # combine all site names into one column
  allSites<-c()
  for(i in 1:numberOfSamples){
    allSites<-c(allSites,sampleList[[i]][,1])
  }
  
  # keep only unique probes
  allSites<-unique(allSites)
  
  # create RRBSmatrix
  RRBSmatrix<-as.data.frame(matrix(ncol=(1+2*numberOfSamples),
                                   nrow=length(allSites)))
  RRBSmatrix[,1]<-allSites
  colnames(RRBSmatrix)[1]<-"chr_start"
  colnames(RRBSmatrix)[2:(numberOfSamples+1)]<-paste(names,"m",sep="_")
  colnames(RRBSmatrix)[(numberOfSamples+2):(2*numberOfSamples+1)]<-paste(names,"cov",sep="_")
  
  # match samples to probes
  for(i in 1:numberOfSamples){
    k<-match(sampleList[[i]][,1],RRBSmatrix[,1])
    RRBSmatrix[k,i+1]<-sampleList[[i]][,2]
    RRBSmatrix[k,i+1+numberOfSamples]<-sampleList[[i]][,3]
  }
  return(RRBSmatrix)
}