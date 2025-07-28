## Getting bumphunter to work

# Dorothea Seiler Vellame 16-11-2018

#BiocInstaller::biocLite("bumphunter")
library(bumphunter)


## load Isabels rTg4510 data
# load filtered data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixFiltered.Rdata")
#RRBSall<-RRBS
#RRBS<-RRBSall[1:100,]
#RRBS<-RRBSall

# remove all NAs (i think this is the problem)
k<-c()
for (i in 1:nrow(RRBS)){
  if (sum(is.na(RRBS[i,2:63]))==0){
    k<-cbind(k,i)
  }
  
}

# subset for testing
RRBS<-RRBS[k,]

# create chr, pos and cl
chrPosList<-strsplit(as.character(RRBS[,1]),"_")
chrPosUnlist<-unlist(chrPosList)

chr<-chrPosUnlist[seq(1,length(chrPosUnlist),2)]
pos<-chrPosUnlist[seq(2,length(chrPosUnlist),2)]

chr1<-chr
pos1<-pos

chr<-chr1[11355:11900]
pos<-pos1[11355:11900]

chr<-chr1
pos<-as.numeric(pos1)


cl<-clusterTest(pos,chr,maxGap=300)


# debug clusterMaker # added line to keep x numeric
clusterTest<-function (chr, pos, assumeSorted = FALSE, maxGap = 300) 
{
  nonaIndex <- which(!is.na(chr) & !is.na(pos))
  Indexes <- split(nonaIndex, chr[nonaIndex]) # seperates by chromosome
  clusterIDs <- rep(NA, length(chr))
  LAST <- 0
  for (i in seq(along = Indexes)) {
    Index <- Indexes[[i]]
    x <- pos[Index]
    if (!assumeSorted) {  # sort the index list
      Index <- Index[order(x)]
      x <- pos[Index]
    }
    y <- as.numeric(diff(x) > maxGap) # x forced to be numeric
    z <- cumsum(c(1, y))
    clusterIDs[Index] <- z + LAST
    LAST <- max(z) + LAST
  }
  clusterIDs
}





# create design matrix
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")

design<-model.matrix(~as.factor(phenotype$Age_months)+phenotype$Genotype+as.factor(phenotype$Age_months)*phenotype$Genotype)


tab<-bumphunter(as.matrix(RRBS[,2:63]),design,chr,pos,cl,cutoff=0.5)

lapply(tab,head)

save(tab,RRBS,file="rTg4510BHRestult.Rdata")

# load tab and turn it into a GRanges object
load("rTg4510BHRestult.Rdata")
library("GenomicRanges")
regions <- GRanges(seqnames = tab$table$chr, 
                                           IRanges(start = tab$table$start, end = tab$table$end),
                                           strand = '*', value = tab$table$value, area = tab$table$area, 
                                           cluster = tab$table$cluster, L = tab$table$L, clusterL = tab$table$clusterL)

## Assign chr lengths
data(hg19Ideogram, package = 'biovizBase')
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]

## Load regionReport
BiocInstaller::biocLite("regionReport")
library(regionReport)
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
## Registered S3 method overwritten by 'dplyr':
##   method               from  
##   as.data.frame.tbl_df tibble
## Make it so that the report will be available as a vignette
original <- readLines(system.file('regionExploration', 'regionExploration.Rmd',
                                  package = 'regionReport'))
vignetteInfo <- c(
  'vignette: >',
  '  %\\VignetteEngine{knitr::rmarkdown}',
  '  %\\VignetteIndexEntry{Basic genomic regions exploration}',
  '  %\\VignetteEncoding{UTF-8}'
)
new <- c(original[1:16], vignetteInfo, original[17:length(original)])
writeLines(new, 'regionReportBumphunter.Rmd')

## Now create the report
report <- renderReport(regions, 'Example bumphunter', pvalueVars = NULL,
                       densityVars = c('Area' = 'area', 'Value' = 'value',
                                       'Cluster Length' = 'clusterL'), significantVar = NULL,
                       output = 'bumphunterExampleOutput', outdir = '.',
                       template = 'regionReportBumphunter.Rmd', device = 'png')