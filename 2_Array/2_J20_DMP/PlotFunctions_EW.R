## Functions for regional plotting 
# Aisha Dahir (A.N.Dahir@exeter.ac.uk)

# Plot DMR Genotypewithin cluster 28/05/2020
plotDMRGenotype <- function(region = region , cols = cols, groups = groups, predictedMeth = predictedMeth, DMR_df = DMR_df){
  
  # This function selects the DMR of interest and plots this including the cluster - the region is highlighted grey
  # DMR = which DMR cluster do you want - you can decide this from the DMR outputs table
  # cols = vector of what colour you want your groups to be in characters
  # groups = groups you want to color should follow the same order as color
  # predictedMeth = output after smoothing methylation
  # DMR_df = Formal class GRanges that can be used to select DMR of interest
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  library(BiSeq)
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]
  
  # Select Cluster where this DMR lies
  #check cluster id
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  #add color to the pheno file for plotting
  
  if (!identical(rownames(pheno), colnames(betas_cluster))) {
    stop("Rownames of betas and pheno matrices are not identical")
  }
  #Add color to pheno file
  pheno$Color <- rep(NA, nrow(pheno))
  #darkslategray1 - TG #black - WT
  pheno$Genotype = as.character(pheno$Genotype)
  pheno$Color = ifelse(pheno$Genotype == groups[1], cols[1], cols[2])
  
  #merge the coordinates and betas for plotting
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  
  ## Plotting
  plot(betas_coord_clus[,"start"], betas_coord_clus[,1],type="o",
       col = "white", ylim = c(0,1),
       xlab = unique(betas_coord_clus$seqnames), ylab = "Smoothed Methylation",
       main = paste( "Cluster: ",unique(betas_coord_clus$cluster.id), sep = ""))
  rect(xleft=DMR$start, xright=DMR$end, ybottom=0.0, ytop=1.0, col="gray70", lty = 0)
  for(i in 1:nrow(pheno)) {
    lines(betas_coord_clus[,"start"],betas_coord_clus[,i], col = pheno[i,"Color"],type="o")
  }
  # par(xpd=TRUE)
  # legend(x ="topright",legend = c("WT", "TG"),col = c(cols[1], cols[2]),
  #        pch = 19, cex = .7,inset=c(0,0))
  add_legend("topright", legend = c(groups[1], groups[2]),col = c(cols[1], cols[2]), pch=20,
             cex=0.8)
}

# Plot DMR Genotype*Age 04/06/2020
plotDMRInteraction <- function(region = region , cols = cols, pch = pch, predictedMeth = predictedMeth, DMR_df = DMR_df){
  
  # This function selects the DMR of interest and plots this including the cluster - the region is highlighted grey
  # DMR = which DMR cluster do you want - you can decide this from the DMR outputs table
  # cols = vector of what colour you want your AGE groups to be in characters
  # pch = vector of line types that your GENOTYPE groups to be in characters
  # predictedMeth = output after smoothing methylation
  # DMR_df = Formal class GRanges that can be used to select DMR of interest
  
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  # Select Cluster where this DMR lies
  #check cluster id
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  #add color to the pheno file for plotting
  
  if (!identical(rownames(pheno), colnames(betas_cluster))) {
    stop("Rownames of betas and pheno matrices are not identical")
  }
  #Add color to pheno file
  pheno$Color <- rep(NA, nrow(pheno))
  #darkslategray1 - TG #black - WT
  # pheno$Color = ifelse(pheno$Genotype == "TG", cols[1], cols[2])
  pheno$Color = ifelse(pheno$Age_months == 2, cols[1],  ifelse(pheno$Age_months == 4, cols[2], ifelse(pheno$Age_months == 6, cols[3], cols[4])))
  
  #Add line types to pheno file
  pheno$LineType <- rep(NA, nrow(pheno))
  # pheno$LineType = ifelse(pheno$Age_months == 2, pch[1],  ifelse(pheno$Age_months == 4, pch[2], ifelse(pheno$Age_months == 6, pch[3], pch[4])))
  pheno$LineType = ifelse(pheno$Genotype == "TG", pch[1], pch[2])
  
  #merge the coordinates and betas for plotting
  
  if (!identical(rownames(betas_cluster), rownames(coordinates_cluster))) {
    stop("Rownames of betas and coordinate matrices are not identical")
  }
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  
  ## Plotting
  plot(betas_coord_clus[,"start"], betas_coord_clus[,1],type="o", pch =  pheno[i,"LineType"],
       col = "white", ylim = c(0,1),
       xlab = unique(betas_coord_clus$seqnames), ylab = "Smoothed Methylation",
       main = paste( "Cluster: ",unique(betas_coord_clus$cluster.id), sep = ""))
  rect(xleft=DMR$start, xright=DMR$end, ybottom=0.0, ytop=1.0, col="gray70", lty = 0)
  for(i in 1:nrow(pheno)) {
    lines(betas_coord_clus[,"start"],betas_coord_clus[,i], col = pheno[i,"Color"],type="o", pch =  pheno[i,"LineType"])
  }
  # par(xpd=TRUE)
  # legend(x ="topright",legend = c("WT_2m", "WT_4m", "WT_6m", "WT_8m",
  #                                 "TG_2m", "TG_4m", "TG_6m", "TG_8m"), 
  #        col = c(cols[1], cols[2], cols[3], cols[4]),
  #        pch = rep(c(1,2), each=  4), cex = .7,inset=c(0,-0.1))
  add_legend("topright", legend = c("WT_2m", "WT_4m", "WT_6m", "WT_8m",
                                    "TG_2m", "TG_4m", "TG_6m", "TG_8m"), 
             col = c(cols[1], cols[2], cols[3], cols[4]),
             pch = rep(c(1,2), each=  4),
             cex=0.8)
  
}

# Plot upstream and downstream to cluster

## Plot the average of the genotype samples within a cluster 08/06/2020
plotDMRGenotypeMedian <- function(region = region , cols = cols, groups = groups, predictedMeth = predictedMeth, DMR_df = DMR_df){
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  library(BiSeq)
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]
  
  # Select Cluster where this DMR lies
  #check cluster id
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  #add color to the pheno file for plotting
  
  if (!identical(rownames(pheno), colnames(betas_cluster))) {
    stop("Rownames of betas and pheno matrices are not identical")
  }
  #Add color to pheno file
  pheno$Color <- rep(NA, nrow(pheno))
  #darkslategray1 - TG #black - WT
  pheno$Genotype = as.character(pheno$Genotype)
  pheno$Color = ifelse(pheno$Genotype == groups[1], cols[1], cols[2])
  
  
  #Pull the WT and TG samples
  WT_ID <- rownames(pheno[which(pheno$Genotype == groups[1]),])
  TG_ID <- rownames(pheno[which(pheno$Genotype == groups[2]),])
  
  betas_WT <- betas_cluster[,WT_ID]
  betas_TG <- betas_cluster[,TG_ID]
  
  median_WT <- apply(betas_WT, 1, median, na.rm=TRUE)
  median_TG <- apply(betas_TG, 1, median, na.rm=TRUE)
  
  #merge the coordinates and betas for plotting
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  
  betas_coord_clus$median_WT <- median_WT
  betas_coord_clus$median_TG <- median_TG
  
  ## Plotting
  par(mar =  c(5.1, 4.1, 4.1, 1))
  plot(betas_coord_clus[,"start"], betas_coord_clus[,1],
       col = "white", ylim = c(0,1),
       xlab = unique(betas_coord_clus$seqnames), ylab = "Smoothed Methylation",
       main = paste( "Cluster: ",unique(betas_coord_clus$cluster.id), sep = ""))
  rect(xleft=DMR$start, xright=DMR$end, ybottom=0.0, ytop=1.0, col="gray70", lty = 0)
  for(i in 1:nrow(pheno)) {
    points(betas_coord_clus[,"start"],betas_coord_clus[,i], col = pheno[i,"Color"])
  }
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_WT"], col = cols[1], lty=2)
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_TG"], col = cols[2], lty=2)
  
  add_legend("topright", legend = c(groups[1], groups[2]),col = c(cols[1], cols[2]), pch=20,
             cex=0.8)
}


## Plot the average of each age in the interaction within a cluster 08/06/2020
plotDMRInteractionMedian <- function(region = region , cols = cols, groups = groups, lty = lty,
                                     pch = pch, predictedMeth = predictedMeth, DMR_df = DMR_df, ages = ages){
  
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  # Select Cluster where this DMR lies
  #check cluster id
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  #add color to the pheno file for plotting
  
  if (!identical(rownames(pheno), colnames(betas_cluster))) {
    stop("Rownames of betas and pheno matrices are not identical")
  }
  #Add color to pheno file
  pheno$Color <- rep(NA, nrow(pheno))
  #darkslategray1 - TG #black - WT
  # pheno$Color = ifelse(pheno$Genotype == "TG", cols[1], cols[2])
  pheno$Color = ifelse(pheno$Age_months == ages[1], cols[1],
                       ifelse(pheno$Age_months == ages[2], cols[2],
                              ifelse(pheno$Age_months == ages[3], cols[3], cols[4])))
  
  #Add line types to pheno file
  pheno$LineType <- rep(NA, nrow(pheno))
  # pheno$LineType = ifelse(pheno$Age_months == 2, pch[1],  ifelse(pheno$Age_months == 4, pch[2], ifelse(pheno$Age_months == 6, pch[3], pch[4])))
  pheno$LineType = ifelse(pheno$Genotype == groups[1], pch[1], pch[2])
  
  # write function to calculate the median of samples of interset
  SamplesMedian <- function(Genotype = Genotype, Age = Age, betas_cluster = betas_cluster){
    samples <- rownames(pheno[which(pheno$Genotype == Genotype & pheno$Age_months == Age),])
    betas_samples <-  betas_cluster[,samples]
    median_samples <- apply(betas_samples, 1, median, na.rm=TRUE)
    
    return(median_samples)
  }
  
  median_WT_2 <- SamplesMedian("WT", ages[1], betas_cluster)
  median_TG_2 <- SamplesMedian("TG", ages[1], betas_cluster)
  median_WT_4 <- SamplesMedian("WT", ages[2], betas_cluster)
  median_TG_4 <- SamplesMedian("TG", ages[2], betas_cluster)
  median_WT_6 <- SamplesMedian("WT", ages[3], betas_cluster)
  median_TG_6 <- SamplesMedian("TG", ages[3], betas_cluster)
  median_WT_8 <- SamplesMedian("WT", ages[4], betas_cluster)
  median_TG_8 <- SamplesMedian("TG", ages[4], betas_cluster)
  
  median_groups <- cbind(median_WT_2,
                         median_TG_2,
                         median_WT_4,
                         median_TG_4,
                         median_WT_6,
                         median_TG_6,
                         median_WT_8,
                         median_TG_8)
  
  #merge the coordinates and betas for plotting
  if (!identical(rownames(betas_cluster), rownames(coordinates_cluster))) {
    stop("Rownames of betas and coordinate matrices are not identical")
  }
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  betas_coord_clus <- cbind(betas_coord_clus, median_groups)
  
  ## Plotting

  par(mar = c(4,4,5,5.5))
  plot(betas_coord_clus[,"start"], betas_coord_clus[,1],type="o", pch =  pheno[i,"LineType"],
       col = "white", ylim = c(0,1),
       xlab = unique(betas_coord_clus$seqnames), ylab = "Smoothed Methylation",
       main = paste( "Cluster: ",unique(betas_coord_clus$cluster.id), sep = ""))
  rect(xleft=DMR$start, xright=DMR$end, ybottom=0.0, ytop=1.0, col="gray70", lty = 0)
  for(i in 1:nrow(pheno)) {
    points(betas_coord_clus[,"start"],betas_coord_clus[,i],col = pheno[i,"Color"], pch =  pheno[i,"LineType"])
  }
  
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_WT_2"], col = "black",  lty=lty[1])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_TG_2"], col = "black",  lty=lty[2])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_WT_4"], col = "brown",  lty=lty[1])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_TG_4"], col = "brown",  lty=lty[2])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_WT_6"], col = "blue",   lty=lty[1])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_TG_6"], col = "blue",   lty=lty[2])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_WT_8"], col = "purple", lty=lty[1])
  lines(betas_coord_clus[,"start"],betas_coord_clus[,"median_TG_8"], col = "purple", lty=lty[2])
  
  
  legend(x ="topright",
          legend = c(paste("WT_", ages[1], sep = ""),
                     paste("WT_", ages[2], sep = ""),
                     paste("WT_", ages[3], sep = ""),
                     paste("WT_", ages[4], sep = ""),
                     paste("TG_", ages[1], sep = ""),
                     paste("TG_", ages[2], sep = ""),
                     paste("TG_", ages[3], sep = ""),
                     paste("WT_", ages[4], sep = "")),
          col = c(cols[1], cols[2], cols[3], cols[4]),
          pch = rep(c(1,2), each=  4),
          lty = rep(c(1,6), each=  4),
         cex = 0.8,
         inset=c(-0.2,-0.2), xpd=TRUE)
  # 
  # add_legend("topright", legend = c("WT_2m", "WT_4m", "WT_6m", "WT_8m",
  #                                   "TG_2m", "TG_4m", "TG_6m", "TG_8m"),
  #            col = c(cols[1], cols[2], cols[3], cols[4]),
  #            pch = rep(c(1,2), each=  4),
  #            lty = rep(c(1,6), each=  4),
  #            cex=0.8)
  
}


# Plot DMR pathology 08/06/2020 - STILL NEEDS CLEANING UP
plotDMRPathology <- function(region = region , cols = cols, groups = groups, lty = lty,
                             pch = pch, predictedMeth = predictedMeth, DMR_df = DMR_df){
  library(BiSeq)
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]
  
  # Select Cluster where this DMR lies
  #check cluster id
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  
  #merge the coordinates and betas for plotting
  
  if (!identical(rownames(betas_cluster), rownames(coordinates_cluster))) {
    stop("Rownames of betas and coordinate matrices are not identical")
  }
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  
  #Add line types to pheno file
  pheno$LineType <- rep(NA, nrow(pheno))
  # pheno$LineType = ifelse(pheno$Age_months == 2, pch[1],  ifelse(pheno$Age_months == 4, pch[2], ifelse(pheno$Age_months == 6, pch[3], pch[4])))
  pheno$LineType = ifelse(pheno$Genotype == groups[1], pch[1], pch[2])
  
  ## Plotting
  ###
  
  path_cols <- colorRampPalette(c(cols[1], cols[2]))(length(which(pheno$ECX >= 0)))
  # par(mar = c(4,4,4,2)) #c(bottom, left, top, right)
  # par(fig=c(0.01,0.90,0.05,0.99)) #c(x1, x2, y1, y2)
  opar <- par(fig=c(0.01,0.90,0.05,0.99),mar = c(4,4,4,2), new=TRUE)
  on.exit(par(opar))
  plot(betas_coord_clus[,"start"], betas_coord_clus[,1],type="o",
       col = "white", ylim = c(0,1),
       xlab = unique(betas_coord_clus$seqnames), ylab = "Smoothed Methylation",
       main = paste( "Cluster: ",unique(betas_coord_clus$cluster.id), sep = ""))
  rect(xleft=DMR$start, xright=DMR$end, ybottom=0.0, ytop=1.0, col="gray70", lty = 0)
  for(i in 1:nrow(pheno)) {
    lines(betas_coord_clus[,"start"],betas_coord_clus[,i], col = path_cols[i],type="o", pch =  pheno[i,"LineType"])
  }
  
  # par(fig=c(0.8,1,0,1), new=TRUE)
  color.bar <- function(lut = lut, min = min, max = max, nticks, title=''){
    ticks=seq(0, (length(which(pheno$ECX >= 0))), len=nticks)
    scale = (length(lut)-1)/(max-min)
    dev.new(width=1.75, height=5)
    plot(c(0,10), c(0,max(pheno$ECX, na.rm = TRUE)), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }
  opar <- par(fig=c(0.8,1,0,1),  new=TRUE)
  on.exit(par(opar))
  color.bar(colorRampPalette(c(cols[1], cols[2]))(round(length(which(pheno$ECX >= 0))/2)),
            min = 00, max = round(length(which(pheno$ECX >= 0))/2), nticks = 11, title = "Pathology")
  
  
  # Add legend for circle and triangle
  par(fig=c(0.85,1,0.05,0.2), new=TRUE)  #c(x1, x2, y1, y2)
  legend(x ="bottomleft",
         legend = c("WT", "TG"),
         pch = c(1,6),
         cex = 1,
         inset=c(-2.5,1.2), 
         xpd=TRUE)
} 


# Plot Boxplots for DMP 18/06/2020
plotDMP <- function(region = 1, position = position, predictedMeth = predictedMeth, 
                    cols = cols, groups = groups,  DMR_df = DMR_df){
  
  library(BiSeq)
  library(ggplot2)
  library(reshape2)
  # region = which DMR cluster do you want - you can decide this from the DMR outputs table
  # position = which position in DMR - can know this by looking at DMR plots or BetaResults for the cluster of interest
  # cols = vector of what colour you want your groups to be in characters
  # groups = groups you want to color should follow the same order as color
  # predictedMeth = output after smoothing methylation
  # DMR_df = Formal class GRanges that can be used to select DMR of interest
  
  
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  DMR <- DMR_df[region,]

  
  Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
  cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
  
  #pull out cluster from coordinates and betas
  coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
  betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
  
  #Add color to pheno file
  pheno$Color <- rep(NA, nrow(pheno))
  pheno$Genotype = as.character(pheno$Genotype)
  pheno$Color = ifelse(pheno$Genotype == groups[1], cols[1], cols[2])
  
  #merge the coordinates and betas for plotting
  betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
  
  #Filter for regions than clusters
  betas_coord_reg <- betas_coord_clus[betas_coord_clus[,"start"] %in% seq(from = DMR_df[region,"start"], to = DMR_df[region,"end"], by = 1),]
  betas_coord_reg2 <- betas_coord_reg[position,c(1:61)] # keep samples and start pos
  betas_coord_reg3 <- cbind(t(betas_coord_reg2), pheno$Genotype, as.numeric(as.character(pheno$Age_months)))
  betas_coord_reg4 <- as.data.frame(betas_coord_reg3)
  colnames(betas_coord_reg4) <- c("Methylation", "Genotype", "Age" )
  betas_coord_reg4$Methylation <- as.numeric(as.character(betas_coord_reg4$Methylation))
  betas_coord_reg4$Genotype <- as.character(betas_coord_reg4$Genotype)
  betas_coord_reg4$Age <- as.numeric(as.character(betas_coord_reg4$Age))
  betas_coord_reg4$Group <- paste(betas_coord_reg4$Genotype, betas_coord_reg4$Age, sep = "_")
  betas_coord_reg4$Group <- factor(betas_coord_reg4$Group, levels = c("WT_2", "TG_2",
                                                                      "WT_4", "TG_4",
                                                                      "WT_6", "TG_6",
                                                                      "WT_8", "TG_8"))
  
  ggplot(betas_coord_reg4, aes(x = Group, y = Methylation)) +
    geom_boxplot(aes(fill = Genotype)) +
    scale_fill_manual(values=c( "#00AEC9", "#000000")) +
    labs(title = paste(betas_coord_reg[1,"seqnames"], betas_coord_reg[1,"start"], sep = "_"),
         x = "Genotype Age Groups",
         y = "Smoothes Methylation") + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Plots for sites of interest
plotDMPGenotype <- function(cpg, color, betaResults, betas, pheno){
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  # if(!identical(pheno$Basename, colnames(betas))){
  #   stop("Pheno basename and Betas colnames are not identical")
  # }
  data_plot <- cbind(coldata[,c("Age_months","Genotype")], t(data[site,]))
  data_plot <- data_plot[complete.cases(data_plot), ]
  par(mar=c(5,5,5,3))
  plot(y=data_plot[,site], x=as.numeric(data_plot$Genotype),
       col = c("black", color)[data_plot$Genotype],
       pch = 19,
       cex=2,
       xlim=c(0.5,2.5),
       ylim=c(min(betas[cpg,]), max(betas[cpg,])),
       axes=FALSE, ann=FALSE)
  axis(side=1, at=1:2, labels=c("WT", "TG"), cex.axis=2)
  axis(2, cex.axis=2)
  title(main = paste(gene, " ", site), ylab = "Methylation (%)", cex.main=2, cex.lab=2)
  
}


plotDMPInteraction <- function(cpg, color, betaResults, betas, pheno){
  
  #cpg = cpg of interest
  #color = color based on model, WT is black by default
  #betaResults = beta results matrix which contains the gene_symbol
  #betas = beta matrix with methylation
  #pheno = phenotype data
  #if(!identical(pheno$Basename, colnames(betas))){
   # stop("Pheno basename and Betas colnames are not identical")
  #}
  
  data_plot <- cbind(pheno[,c("Age_months","Genotype")], t(betas[cpg,]))
  data_plot$col <-  ifelse(data_plot$Genotype == "WT","#000000", color)
  ages <- sort(unique(pheno$Age_months))
  medians_WT <- c()
  medians_TG <- c()
  for (age in ages) {
    medians_WT <- c(medians_WT, median(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="WT"))[,cpg]))
    medians_TG <- c(medians_TG, median(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="TG"))[,cpg]))
  }
  
  par(mar=c(5,5,5,3))
  plot(y=betas[cpg,], x=as.numeric(pheno$Age_months),
       col = data_plot$col,
       pch = 19,
       cex=2,
       xlim=c(min(ages)-0.5,max(ages)+0.5), #make the x axis range depend on the age interval
       ylim=c(min(betas[cpg,]), max(betas[cpg,])),
       axes=FALSE, ann=FALSE)
  axis(side=1, at=ages, labels=as.character(ages), cex.axis=2)
  axis(2, cex.axis=2)
  title(main = paste(cpg, betaResults[cpg, "Gene_Symbol"], sep = " "), 
        xlab = "Age (months)", ylab = "Methylation (%)", cex.main=2, cex.lab=2)
  lines(ages, medians_WT, lty=2, col="black", lwd=3)
  lines(ages, medians_TG, lty=2, col= color, lwd=3)
  
}


plotDMPPathology <- function(cpg, color, betaResults, betas, pheno){
  
  if(!identical(pheno$Basename, colnames(betas))){
    stop("Pheno basename and Betas colnames are not identical")
  }
  
  coldata_pathology <- pheno
  data <- betas[, coldata_pathology$Basename]
  coldata <- coldata_pathology
  identical(coldata_pathology$Basename, colnames(data))
  ages <- sort(unique(pheno$Age_months))
  index <- cpg
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  data_plot <- cbind(coldata[,c("Age_months","Genotype")], t(data[site,]))
  data_plot <- data_plot[complete.cases(data_plot), ]
  medians_WT <- c()
  medians_TG <- c()
  for (age in unique(coldata$Age_months)) {
    medians_WT <- c(medians_WT, median(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="WT"))[,site]))
    medians_TG <- c(medians_TG, median(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="TG"))[,site]))
  }
  
  par(mar=c(5,5,5,3))
  plot(y=data_plot[,site], x=as.numeric(data_plot$Age_months),
       col = c("black", color)[data_plot$Genotype],
       pch = 19,
       cex=2,
       xlim=c(min(ages)-0.5,max(ages)+0.5), #make the x axis range depend on the age interval
       ylim=c(min(data_plot[,site]), max(data_plot[,site])),
       axes=FALSE, ann=FALSE)
  axis(side=1, at=ages, labels=as.character(ages), cex.axis=2)
  axis(2, cex.axis=2)
  title(paste(gene, " ", site), xlab = "Age (months)", ylab = "Methylation (%)", cex.main=2, cex.lab=2)
  lines(ages, medians_WT, lty=2, col="black", lwd=3)
  lines(ages, medians_TG, lty=2, col=color, lwd=3)
}


### plots genotypes across two tissues

plotTissueGenotypeDMP <- function(cpg, betaResults, betas, pheno, colour, column){
  # column is the term of interset or the column with the fdr p value you want to write
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue")], data[site,])
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  p <- ggplot(data_plot, aes(x=Genotype, y=data_plot[,4], color=Genotype)) + 
    geom_point() +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", colour)) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - "),
         subtitle =paste("FDR Pvalue: ", signif(betaResults[i, column],3), sep = ""))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}


plotPathologyDMP <- function(cpg, betaResults, betas, pheno, colour, column){
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue","Pathology")], t(data[site,]))
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue","Pathology", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  
  p <- ggplot(data_plot, aes(x=Pathology, y=data_plot[,5], color=Genotype)) + 
    geom_point() +
    #facet_wrap(~Tissue, scales = "free") +
    scale_color_manual(values=c("black", colour)) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - "),
         subtitle =paste("FDR Pvalue: ", signif(betaResults[i, column],3), sep = "")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

  print(p)
}



# Plots for sites of interest
plotDMPAge <- function(cpg, color, betaResults, betas, pheno){
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  # if(!identical(pheno$Basename, colnames(betas))){
  #   stop("Pheno basename and Betas colnames are not identical")
  # }
  data_plot <- cbind(coldata[,c("Age_months","Genotype")], t(data[site,]))
  data_plot <- data_plot[complete.cases(data_plot), ]
  par(mar=c(5,5,5,3))
  plot(y=data_plot[,site], x=as.numeric(data_plot$Age_months),
       #col = c("black", color)[data_plot$Age_months],
       pch = 19,
       cex=2,
       #xlim=c(0.5,2.5),
       ylim=c(min(betas[cpg,]), max(betas[cpg,])), ann = FALSE)
  #axes=FALSE, ann=FALSE)
  #axis(side=1, at=1:2) #labels=c("WT", "TG"), cex.axis=2)
  #axis(2, cex.axis=2)
  title(main = paste(gene, " ", site), ylab = "Methylation (%)", cex.main=2, cex.lab=2, xlab = "Age (months)")
  
}
