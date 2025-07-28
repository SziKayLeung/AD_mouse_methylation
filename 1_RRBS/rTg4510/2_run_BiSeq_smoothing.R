#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## run biseq from bismark coverage output
## run biseq smoothing
## 11.06.2024: rTg4510 all samples
## --------------------------------

message("Using cores: ", parallel::detectCores())

## ---------- packages -----------------

#BiocManager::install("BiSeq")
suppressMessages(library("BiSeq"))
suppressMessages(library("stringr"))
suppressMessages(library("parallel"))
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/import.config")

model = "rTg4510"
## --------------- Import bismark output -------------------- 

# list <bismark.cov.gz>
message("Import bismark coverage files after merging CpG")
infile <- list.files(dirnames$rTg4510bismark, pattern = "merged_CpG_evidence.cov", all.files=T, full.names=T)
print(infile)

# subset to samples with coverage
samples <- word(unlist(lapply(infile, function(x) basename(x))),c(2), sep = fixed("_"))
if(length(setdiff(row.names(phenotype$rTg4510), samples)) == 0){
  print("Processing all samples")
}else{
  message("Missing coverage output for following samples:")
  setdiff(row.names(phenotype$rTg4510), samples)
}

# import
rrbs <- readBismark(files = infile, colData = phenotype$rTg4510)
message("Bismark files fully imported")


## --------------- Definition of CpG clusters -------------------- 

# # Within a BSraw object clusterSites searches for agglomerations of CpG sites across all samples.
# # In a first step the data is reduced to CpG sites covered in round(perc.samples*ncol(object)) samples (here:  5 samples), these are called ’frequently covered CpG sites’.
# # In a second step regions are detected where not less than min.sites frequently covered CpG sites are sufficiantly close to each other (max.dist).

rrbs.clust.unlim <- clusterSites(object = rrbs,
                                 groups = colData(rrbs)$Genotype,
                                 perc.samples = 0.5,
                                 min.sites = 5,
                                 max.dist = 500)
# max.dist: check if output makes sense and perhaps change depending on distance of CpGs in output ##############################################
# rrbs.clust.unlimis a againBSrawobject but restricted to CpG sites withinCpG clusters.  Each CpG site is assigned to a cluster.

# # The underlying CpG clusters can also be converted to aGRangesobject withthe start and end positions
#GRangesobject <- clusterSitesToGR(rrbs.clust.unlim)

message("created rrbs.clust.unlim and GRanges")

## --------------- Smoothing methylation data -------------------- 

# In the smoothing step CpG sites with high coverages get high weights. To reduce bias due to unusually high coverages we limit the coverage, e.g. to the 90% quantile.
ind.cov <- totalReads(rrbs.clust.unlim) > 0
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)

# We then smooth the methylation values of CpG sites within the clusters with the default bandwidthh = 80 base pairs.
# It is possible - and recommended- to parallelize this step by settingmc.cores, to 6 cores for instance, if there are 6 available.
#message("Predicting methylation")
#message("Using cores: ", parallel::detectCores())
#predictedMeth <- predictMeth(object = rrbs.clust.lim, mc.cores = parallel::detectCores())
## predictedMeth is a BSrel object with smoothed relative methylation levels for each CpG site within CpG clusters

# Assuming bsraw is your BSraw object
total_rows <- nrow(rrbs.clust.lim)
chunk_size <- 100
num_chunks <- ceiling(total_rows / chunk_size)

predictMethChunk <- function(start_row, end_row) {
  chunk_data <- rrbs.clust.lim[start_row:end_row, ]
  predictMeth(object = chunk_data, mc.cores = 6)  # Use 1 core per chunk to avoid overloading
}

# Initialize list to store results
predictedMeth_list <- vector("list", num_chunks)

# Function to process each chunk and return data frame
processChunk <- function(i) {
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, total_rows)
  
  cat(sprintf("Processing chunk %d/%d (%d to %d)\n", i, num_chunks, start_row, end_row))
  
  # Predict methylation for current chunk
  predictedMeth_chunk <- predictMethChunk(start_row, end_row)
  
  # Convert to data frame and assign row names
  meth_df <- as.data.frame(methLevel(predictedMeth_chunk))
  coordinates <- as.data.frame(predictedMeth_chunk@rowRanges)
  rownames(meth_df) <- paste0(coordinates$seqnames, ":", coordinates$start)
  
  cat(sprintf("Finished chunk %d/%d\n", i, num_chunks))
  
  return(meth_df)
}

# Process chunks in parallel and combine results
predictedMeth_list <- lapply(1:num_chunks, processChunk)
RRBS_smoothbetas <- do.call(rbind, predictedMeth_list)

# save output
save(predictedMeth_list, file = paste0(dirnames$processed,"/rrbs/",model,"_predictedMeth_list.RData"))
save(RRBS_smoothbetas, file = paste0(dirnames$processed,"/rrbs/",model,"_RRBS_SmoothBetas.RData"))
save(rrbs.clust.lim, file = paste0(dirnames$processed,"/rrbs/biseq_rrbs_clust_lim_",model,".RData"))
save(rrbs, file = paste0(dirnames$processed,"/rrbs/biseq_raw_",model,".RData"))

## --------------- Merge smoothed DNA methylation (Biseq) and raw methylation, so we don't lose sites that were not asigned to clusters -------------------- 

# Filter raw methylation data by coverage
#dim(rrbs)[1] # 3,663,442 sites 
#rrbs.reduced <- filterBySharedRegions(object=rrbs, groups=colData(rrbs)$Genotype, perc.samples=0.5, minCov=10)

# convert to methylation values
#rrbs.rel <- rawToRel(rrbs.reduced) 
#RRBS_rawbetas <- as.data.frame(methLevel(rrbs.rel))
#nrow(RRBS_rawbetas) # 1,125,882
#coordinates <- as.data.frame(rrbs.reduced@rowRanges)
#rownames(RRBS_rawbetas) <- paste0(coordinates$seqnames, ":", coordinates$start)

# Get smooth methylation matrix
#RRBS_smoothbetas <- as.data.frame(methLevel(predictedMeth))
#nrow(RRBS_smoothbetas) # 1,294,995
#coordinates <- as.data.frame(predictedMeth@rowRanges)
#rownames(RRBS_smoothbetas) <- paste0(coordinates$seqnames, ":", coordinates$start)

message("Save output")
# save biseq raw methylation matrix
#save(rrbs, file = paste0(dirnames$processed,"/rrbs/biseq_raw_rTg4510_parallel.RData"))
#save(rrbs.clust.unlim, GRangesobject, predictedMeth, file = paste0(dirnames$processed,"/rrbs/biseq_rTg4510_parallel.RData"))
#save(RRBS_smoothbetas, file = "/lustre/projects/Research_Project-MRC148213/lsl693/rrbs_ad_mice/0_ZenOutput/1_processed/rrbs/rTg4510_RRBS_SmoothBetas_parallel.RData")
