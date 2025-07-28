#Emma Walker
# E.M.Walker@exeter.ac.uk
# 04/02/19


###GRanges for original rdata files

# load packages
library(stringr)
library(GenomicRanges)
library(dplyr)


# get csv files for annotation and format
temp <- list.files(path = "/mnt/data1/Thea/IsabelSamples/data/", pattern = "GenotypeStats.Rdata")



for(i in 1:length(temp)){
  
  # Get file
  inpath <- paste0("/mnt/data1/Thea/IsabelSamples/data/", temp[i])
  load(inpath)
  File <- as.data.frame(RRBSGenotypeStat)
  colnames(File) <- c("F_value", "Degree_of_freedom", "Effect_size", "p_value", "FDR", "n_Total")
  File$X <- row.names(File)
  
  
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
  ref_outpath <- paste0("/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/GRanges_Output/",
                        temp[i], "ref.csv")
  
  write.csv(seg.df, ref_outpath)
  
  
  ########### create final formated file for Isabel
  
  # get orignal csv files and format location
  
  orig <- as.data.frame(RRBSGenotypeStat)
  colnames(orig) <- c("F_value", "Degree_of_freedom", "Effect_size", "p_value", "FDR", "n_Total")
  orig$Location <- row.names(File)
  orig$Location <- gsub("_", ":", orig$Location)
  
  # combine using cbind
  final <- cbind.data.frame(orig, seg.df$Gene, seg.df$Nearest_gene, seg.df$`near.df$distance`)
  colnames(final)[8:10] <- c("Gene", "Nearest_Gene", "Distance_to_gene")
  
  # write final output to csv
  final_outpath <- paste0("/mnt/data1/EmmaW/Isabel/annotate_rrbs_locations/For_Isabel/table_version",
                          temp[i], "annotated.csv")
  
  write.table(final, final_outpath, sep = "\t")
  
}
