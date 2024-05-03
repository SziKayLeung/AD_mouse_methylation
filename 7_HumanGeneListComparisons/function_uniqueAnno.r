
#define unique anno function
uniqueAnno<-function(row){
  if(row != ""){
    return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
  } else {
    return(row)
  }
}