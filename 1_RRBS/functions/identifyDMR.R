identify_DMR <- function(betaResults, vario.sm.null){

  # numeric column for position
  betaResults <- as.data.frame(betaResults)
  betaResults$pos <- as.integer(as.character(betaResults$pos))
  
  # estimate variogram of methylation between two CpG sites within a cluster
  # calculate the z.score
  vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
  vario.sm.null$pValsList <- vario.aux$pValsList

  ##Correlaton is based on z-scores
  locCor <- estLocCor(vario)
  
  # test each CpG cluster for the presence of at least one differentially methylated  location  
  # at q what  can  be  interpreted  as  the  size-weighted  FDR  on clusters:
  clusters.rej <- testClusters(locCor, FDR.cluster = 0.05)
  
  # trim the rejected CpG clusters that is to remove the not dierentially methylated CpG sites 
  # at q1 what can be interpreted as the location-wise FDR
  clusters.trimmed <- trimClusters(clusters.rej, FDR.loc = 0.05)
  
  # distance 500bp
  DMR <- findDMRs(clusters.trimmed, max.dist = 500, diff.dir = TRUE)
  
  return(DMR)
}
