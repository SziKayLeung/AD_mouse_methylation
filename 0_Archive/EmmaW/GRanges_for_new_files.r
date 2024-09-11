#Emma Walker
# E.M.Walker@exeter.ac.uk
# 04/02/19


###GRanges for edited data files

# load packages
library(stringr)
library(GenomicRanges)
library(dplyr)


#for new files 21/02/19
temp <- list.files(path = "/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/", pattern = "NEW.csv", 
                   full.names = T)



for(i in 1:length(temp)){
  
  # Get file
  File <- read.csv(temp[i])

  
  # add columns to csv files for chromosome, location, min location and max location
  File$X <- gsub("M", "MT", File$X)
  File$loc <- sub(".*_", "", File$X)
  File$X <- gsub("_", ":", File$X)
  File$query <-paste0(File$X, "-", File$loc)
  
  # turn File$query into a GRanges object 
  seg <- GRanges(File$query)
  
  # get predownloaded biomart file for reference genome
  mart <- openxlsx::read.xlsx("/mnt/data1/EmmaW/Isabel/mart_export_reformatted.xlsx")
  
  # format turn reference file and turn into GRanges object
  ref <- data.frame()
  ref <-paste0("chr", mart$Chromosome.Name, ":", mart$`Gene.Start.(bp)`, "-", mart$`Gene.End.(bp)`)
  ref <- GRanges(ref)
  mcols(ref)$Genes <- mart$Associated.Gene.Name
  
  # find overlap in query and ref regions
  olap <- findOverlaps(seg, ref)
  
  # add gene annotations to query file
  as(olap, "List")
  mcols(seg)$Gene <- extractList(mcols(ref)$Genes,
                                 as(olap, "List"))
  
  # find name of nearest gene and distance, and add to query file
  near <- distanceToNearest(seg, ref)
  as(near, "List")
  mcols(seg)$Nearest_gene <- extractList(mcols(ref)$Genes, as(near, "List"))
  
  seg.df <- as.data.frame(seg)
  near.df <- as.data.frame(near)
  
  seg.df <- cbind.data.frame(seg.df, near.df$distance)
  
  
  #write annoated file to csv
  
  temp2 <- list.files(path = "/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/", pattern = "NEW.csv")
  ref_outpath <- paste0("/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/GRanges_Output/",
                        temp2[i], "ref.csv")

  
  write.csv(seg.df, ref_outpath)
  
  
  ########### create final formated file for Isabel
  
  # get orignal csv files and format location
  
  orig <- read.csv(temp[i])
  orig$X <- gsub("_", ":", orig$X)
  colnames(orig)[1] <- "Location"
  
  # combine using cbind
  final <- cbind.data.frame(orig, seg.df$Gene, seg.df$Nearest_gene, seg.df$`near.df$distance`)
  colnames(final)[8:10] <- c("Gene", "Nearest_Gene", "Distance_to_gene")
  
  # write final output to csv
  final_outpath <- paste0("/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/For_Isabel/table_version",
                          temp2[i], "annotated.csv")
  
  write.table(final, final_outpath, sep = "\t")
  
}
