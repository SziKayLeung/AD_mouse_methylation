
#functions to create gene list table for report

makeGeneList <- function(overlap_list){
  
  names(overlap_list) <- tab$Model.Term
  
  n.obs <- sapply(overlap_list, length)
  seq.max <- seq_len(max(n.obs))
  
  mat <- as.data.frame(sapply(overlap_list, "[", i = seq.max))
  mat[] <- lapply(mat, as.character)
  mat[is.na(mat)] <- ""
  
  return(mat)
  
}

ECX.makeGeneList <- function(overlap_list){
  
  names(overlap_list) <- tab$Model.Term
  
  n.obs <- sapply(overlap_list, length)
  seq.max <- seq_len(max(n.obs))
  
  mat <- as.data.frame(t(sapply(overlap_list, "[", i = seq.max)))
  mat[] <- lapply(mat, as.character)
  mat[is.na(mat)] <- ""
  
  return(mat)
  
}


