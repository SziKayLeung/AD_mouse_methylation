# Function to compare lists of genes and report on their similarity
## Input: txt files containing a single list of genes each
# i.castanh@bidmc.harvard.edu

# Load packages
library(ggplot2)
library(rcartocolor)
library(reshape2)
library(wesanderson)
library(pheatmap)
library(dplyr)


########################################################################
# Preparation of the files
########################################################################

# Store files containing list of genes to compare in the same directory
## Format: .txt
## Column names: no
## Column 1 = genes to compare (gene symbol)
## File name = what the list represents
### E.g.: EC_MultiregionCIRCUITS_AD-CTRL_Ast.txt
#### Brain region, data set / origin, comparison / facet, cell type


########################################################################
# Function to compare genes
########################################################################

compGenes <- function(directory) {
  
  ## directory = path to directory where lists of genes to compare are stored
  
  
  if(class(directory) == "character") {
    
    
    ################### Prepare the data ################### 
    # Read files with .csv extension in directory of interest
    files_paths <- list.files(path = directory,
                        pattern = "\\.txt$",
                        full.names = TRUE)
    ## Future: read sub-directories (recursive = TRUE)
    
    # Save name of files
    files_names <- list.files(path = directory,
                        pattern = "\\.txt$",
                        full.names = FALSE)
    ## Future: read-sub-directories
    ## (use recursive = TRUE, but include.dirs = FALSE is not working and the directory is also saved in the name of the file)
    ## Tried dir_ls() with no success [library(fs)]
    
    # Save what each file represents (from name of file, that should represent list)
    lists_names <-   sapply(strsplit(files_names, split='.txt', fixed=TRUE),
                            function(x) (x[1]))

    
    ################### Make list of lists ###################
    # Read list of genes in each file in the directory of interest (Future: and sub-directories)
    # and save as a list of lists
    
    ## Create empty list
    my_list <- vector(mode='list', length = length(files_paths))
    
    ## Save genes in list of lists
    for (i in 1:length(files_paths)){
      path_to_file <- files_paths[i]
      temp_df <- read.table(file = path_to_file, sep = "\t")
      temp_names <- unique(temp_df[[1]]) # unique just in case there are genes that are repeated
      my_list[[i]] <- temp_names
      names(my_list)[i] <- lists_names[i]
    }
  }
  
  ################### Calculate the similarity coefficient ###################
  # Compare each element (in the list of lists) with all other elements (and with itself)
  # by calculating the Jaccard index (similarity coefficient)
  # and save as matrix
  
  ## Pull the number of intersecting genes in each comparison
  nGenes_mat <- sapply(my_list,
                          function(x) sapply(my_list,
                                             function(y) length(intersect(x,y))
                                             )
                       )
  
  ### Save matrix containing the number of intersecting genes in a file
  write.csv(nGenes_mat, file = paste(directory,
                                     "/similarityMatrix_numberOfIntersectingGenes.csv", sep = ""))
  
  ### Reshape the matrix
  nGenes_melt <- reshape2::melt(nGenes_mat)
  
  ## Calculate the Jaccard index (similarity coefficient) and save as matrix and data frame
  
  JacInd_mat <- sapply(my_list,
                       function(x) sapply(my_list,
                                          function(y) length(intersect(x,y)) / length(union(x,y))
                                          )
                       )
  
  ### Save matrix containing the Jaccard index for each comparison in a file
  write.csv(JacInd_mat, file = paste(directory,
                                     "/similarityMatrix_JaccardIndex.csv", sep = ""))
  
  ### Reshape the matrix
  JacInd_melt <- reshape2::melt(JacInd_mat)
  
  
  ################### Save list of common genes for each comparison ###################
  # Export the common genes to a file
  
  ### Create empty list (with length = number of comparisons)
  nComparisons <- nrow(nGenes_melt)
  list_export <- vector(mode = "list", length = nComparisons)
  ### Start comparison number from zero
  comparison <- 0
  ### For each element of the list, compare the intersection with all other elements of the list
  for (j in 1:length(my_list)) {
    list1 <- my_list[[j]] # Get the element of the list to compare with all other elements
    for (k in 1:length(my_list)) {
      list2 <- my_list[[k]] # Get each of the other elements to compare
      genes_to_save <- intersect(list1, list2)
      comparison <- comparison + 1
      list_export[[comparison]] <- as.factor(genes_to_save)
      names(list_export)[comparison] <- paste(names(my_list)[j], "_&_", names(my_list)[k], sep = "")
    }
  }
  
  ### Save to file (as .csv to be in a different format compared to the lists to test)
  options(max.print=1000000) # To increase the max number of entries I can print in R to save all genes
  capture.output(list_export, file = paste(directory,
                                           "/ListOverlapGenes_fromSimilarityMatrix.csv", sep = ""))

  
  ################### Plot the results ###################
  # Prepare the data to plot as a diagonal correlation matrix (and not a redundant matrix)
  
  cormat <- JacInd_mat

  ## Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  ## Get upper triangle of matrix
  upper_tri <- get_upper_tri(cormat)
  
  ## Melt the correlation matrix
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
  
  
  # Plot the matrix as a diagonal correlation heatmap: simple
  ggheatmap <- ggplot(data = melted_cormat, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) +
    # geom_text(aes(label = value), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = "dark orange",
                        limit = c(0,1), name = "Jaccard index") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          axis.text = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) 
  
  ggsave(filename = "heatmap1_simple.png", plot = ggheatmap, path = directory, width = 13, height = 13, dpi = 300)
  
  
  # Plot the matrix as a diagonal correlation heatmap: add correlation coefficients
  
  ## color & number = Jaccard Index
  
  ggheatmap2 <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 2.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  ggsave(filename = "heatmap2_JaccardIndex.png", plot = ggheatmap2, path = directory, width = 13, height = 13, dpi = 300)

  
  # Plot the matrix as a diagonal correlation heatmap: add number of genes
  
  ## color = Jaccard Index
  ## number = common genes
  
  # Prepare the matrix showing the number of genes to plot as a diagonal correlation matrix
  
  nGenes_cormat <- nGenes_mat
  
  ## Get upper triangle of matrix
  nGenes_upper_tri <- get_upper_tri(nGenes_cormat)
  
  ## Melt the correlation matrix
  melted_nGenes_cormat <- reshape2::melt(nGenes_upper_tri, na.rm = TRUE)
  
  ggheatmap3 <- ggheatmap + 
    geom_text(data = melted_nGenes_cormat, aes(Var2, Var1, label = value), color = "black", size = 3.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  ggsave(filename = "heatmap3_nGenes.png", plot = ggheatmap3, path = directory, width = 13, height = 13, dpi = 300)
  
  # Adapted using pheatmap to add annotation bars
  annotateGroups <- function(name){
    if(grepl("human", name)){
      return("human")
    }else if(grepl("rTg4510",name)){
      return("rTg4510")
    }else{
      return("J20")
    }
  }
  
  annotateMethod <- function(name){
    if(grepl("array", name)){
      return("Array")
    }else if(grepl("rrbs", name)){
      return("RRBS")
    }else{
      return(NA)
    }
  }
  
  annotateModel <- function(name){
    if(grepl("genotype", name)){
      return("Genotype")
    }else if(grepl("interaction", name)){
      return("Interaction")
    }else if(grepl("pathology", name)){
      return("Pathology")
    }else{
      return(NA)
    }
  }
  
  color_Tg4510_TG <- "#00AEC9"
  color_J20_TG <- "#FF5A62"
  
  my_colour = list(
    species = c(human = "black", rTg4510 = wes_palette("Darjeeling2")[2], J20 = wes_palette("Darjeeling1")[1]),
    approach = c(Array = alpha(wes_palette("Chevalier1")[1],0.6), RRBS = wes_palette("Chevalier1")[1]),
    model = c(Genotype = wes_palette("GrandBudapest2")[1], Interaction = wes_palette("GrandBudapest1")[2], Pathology = wes_palette("GrandBudapest1")[3])
  )
  
  # colour using jaccard index (melted_cormat$Var1), but box the number of common genes
  # create annotation columns of species, approach and model based on file names
  annotation <- data.frame(files = unique(melted_cormat$Var1))
  annotation$species <- sapply(annotation$files, annotateGroups) 
  annotation$approach <- sapply(annotation$files, annotateMethod) 
  annotation$model <- sapply(annotation$files, annotateModel) 
  # replace NA with "" 
  nGenes_upper_tri[is.na(nGenes_upper_tri)] <- ""
  
  # set up pheatmap with rownames
  rownames(annotation) <- annotation$files
  annotation <- annotation %>% select(-files)
  log10upper_tri <- log10(upper_tri+1)
  ggheatmap4 <- pheatmap(upper_tri, cluster_cols = FALSE, cluster_rows = FALSE,annotation_col = annotation,
                        annotation_row = annotation, display_numbers = nGenes_upper_tri, fontsize = 14,
                        annotation_colors = my_colour, color = rev(hcl.colors(30, "BluYl")), border_color = "white",
                        na_col = "white",cellheight = 40, cellwidth = 40, legend = FALSE)
  
  return(ggheatmap4)

}
