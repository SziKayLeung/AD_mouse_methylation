#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS individual plots

# R 3.4.3

setwd("/mnt/data1/isabel/RRBS/")

library(VennDiagram)
#library("DESeq2")
#library("pheatmap")
#library("WGCNA")
#library(vioplot)

color_Tg4510_TG <- "#00AEC9"
color_J20_TG <- "#FF5A62"

# Get full data for rT4510 (all sites, sig or not!)
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeStats.Rdata")
Tg4510_full_data <- as.data.frame(RRBSGenotypeStat)
rownames(Tg4510_full_data) <- sapply(rownames(Tg4510_full_data), gsub, pattern="_",replacement=":")

# Get full data for J20 (all sites, sig or not!)
load(file="/mnt/data1/Thea/IsabelSamples/data/J20GenotypeStats.Rdata")
J20_full_data <- as.data.frame(RRBSGenotypeStat)
rownames(J20_full_data) <- sapply(rownames(J20_full_data), gsub, pattern="_",replacement=":")

# make sure I have the same sites in both datasets
common_sites_fulldata <- intersect(rownames(Tg4510_full_data), rownames(J20_full_data))
Tg4510_data <- Tg4510_full_data[common_sites_fulldata,]
J20_data <- J20_full_data[common_sites_fulldata,]

## SITES
# Significant (FDR <0.05) sites
Tg4510_FDR <- Tg4510_data[which(Tg4510_data[,"FDR_p"] < 0.05),]
J20_FDR <- J20_data[which(J20_data[,"FDR_p"] < 0.05),]

Tg4510_sites <- rownames(Tg4510_FDR)
J20_sites <- rownames(J20_FDR)
common_sites <- intersect(Tg4510_sites, J20_sites)

# Venn diagram #Tg4510_sites #J20_sites
overlap <- calculate.overlap(x = list('rTg4510' = Tg4510_sites, 'J20' = J20_sites))
area_overlap <- sapply(overlap, length)

g <- draw.pairwise.venn(area1 = area_overlap[1],
                        area2 = area_overlap[2],
                        cross.area = area_overlap[3],
                        category = c("rTg4510", "J20"),
                        ind = FALSE,
                        inverted = FALSE,
                        scaled = TRUE,
                        ext.text = TRUE,
                        lwd = 1, 
                        ext.line.lwd = 1,
                        ext.dist = -0.15,
                        ext.length = 0.9,
                        ext.pos = -4,
                        fill = c(color_Tg4510_TG, color_J20_TG),
                        cex = 6,
                        cat.cex= 7,
                        cat.col = c("black", "black"),
                        cat.dist = rep(-0.05, 1), # cat.dist = rep(0.025, 2) # c(-0.16, -0.09)
                        cat.pos = c(5,5),
                        rotation.degree = 35)

# check for switch in labels
if(area_overlap[1] != (as.integer(g[[5]]$label) + as.integer(g[[7]]$label)) && area_overlap[2] !=  (as.integer(g[[6]]$label) + as.integer(g[[7]]$label))){
  # change inverted to TRUE
  g <- draw.pairwise.venn(area1 = area_overlap[1],
                          area2 = area_overlap[2],
                          cross.area = area_overlap[3],
                          category = c("rTg4510", "J20"),
                          ind = FALSE,
                          inverted = TRUE,
                          scaled = TRUE,
                          ext.text = TRUE,
                          lwd = 1, 
                          ext.line.lwd = 1,
                          ext.dist = -0.15,
                          ext.length = 0.9,
                          ext.pos = -4,
                          fill = c(color_Tg4510_TG, color_J20_TG),
                          cex = 6,
                          cat.cex = 7,
                          cat.col = c("black", "black"),
                          cat.dist = rep(-0.05, 1), # cat.dist = rep(0.025, 2) # c(-0.16, -0.09)
                          cat.pos = c(0,5),
                          rotation.degree = 35)
    # switch labels
  tmp_var      <- g[[6]]$label
  g[[6]]$label <- g[[5]]$label
  g[[5]]$label <- tmp_var
  rm(tmp_var)
}

jpeg("Tg4510_vs_J20_sites.jpg", width = 1200, height = 1200)
plot.new()
title(sub = "Total sites = 867012", cex.sub = 6, line = -4, outer = TRUE)
grid.draw(g)
dev.off()


## GENES
Tg4510_FDR_annotations <- read.csv("/mnt/data1/isabel/RRBS/rTg4510GenotypeStats_cleanNEW_annotated.csv", row.names=1, stringsAsFactors=FALSE)
J20_FDR_annotations <- read.csv("/mnt/data1/isabel/RRBS/J20GenotypeStatsNEW_annotated.csv", row.names=1, stringsAsFactors=FALSE)

Tg4510_genes <- unique(Tg4510_FDR_annotations$Nearest_Gene)
J20_genes <- unique(J20_FDR_annotations$Nearest_Gene)
common_genes <- intersect(Tg4510_genes, J20_genes)

# Venn diagram #Tg4510_genes #J20_genes
overlap <- calculate.overlap(x = list('rTg4510' = Tg4510_genes, 'J20' = J20_genes))
area_overlap <- sapply(overlap, length)

g <- draw.pairwise.venn(area1 = area_overlap[1],
                        area2 = area_overlap[2],
                        cross.area = area_overlap[3],
                        category = c("rTg4510", "J20"),
                        ind = FALSE,
                        inverted = FALSE,
                        scaled = TRUE,
                        ext.text = TRUE,
                        lwd = 1, 
                        ext.line.lwd = 1,
                        ext.dist = -0.15,
                        ext.length = 0.9,
                        ext.pos = -4,
                        fill = c(color_Tg4510_TG, color_J20_TG),
                        cex = 6,
                        cat.cex= 7,
                        cat.col = c("black", "black"),
                        cat.dist = rep(-0.025, 2), # cat.dist = rep(0.025, 2) # c(-0.16, -0.09)
                        #cat.pos = c(5,5),
                        rotation.degree = 35)

# check for switch in labels
if(area_overlap[1] != (as.integer(g[[5]]$label) + as.integer(g[[7]]$label)) && area_overlap[2] !=  (as.integer(g[[6]]$label) + as.integer(g[[7]]$label))){
  # change inverted to TRUE
  g <- draw.pairwise.venn(area1 = area_overlap[1],
                          area2 = area_overlap[2],
                          cross.area = area_overlap[3],
                          category = c("rTg4510", "J20"),
                          ind = FALSE,
                          inverted = TRUE,
                          scaled = TRUE,
                          ext.text = TRUE,
                          lwd = 1, 
                          ext.line.lwd = 1,
                          ext.dist = -0.15,
                          ext.length = 0.9,
                          ext.pos = -4,
                          fill = c(color_Tg4510_TG, color_J20_TG),
                          cex = 6,
                          cat.cex = 7,
                          cat.col = c("black", "black"),
                          cat.dist = rep(-0.025, 2), # cat.dist = rep(0.025, 2) # c(-0.16, -0.09)
                          #cat.pos = c(0,5),
                          rotation.degree = 35)
    # switch labels
  tmp_var      <- g[[6]]$label
  g[[6]]$label <- g[[5]]$label
  g[[5]]$label <- tmp_var
  rm(tmp_var)
}

jpeg("Tg4510_vs_J20_genes.jpg", width = 1200, height = 1200)
plot.new()
#title(sub = "Total genes", cex.sub = 6, line = -4, outer = TRUE)
grid.draw(g)
dev.off()

# Save stats for common sites
common_sites_statsTg4510 <- Tg4510_data[common_sites,]
new_common_sites_statsTg4510 <- common_sites_statsTg4510[ order(row.names(common_sites_statsTg4510)), ]
new_Tg4510_FDR_annotations <- Tg4510_FDR_annotations[ order(row.names(Tg4510_FDR_annotations)), ]
new_Tg4510_FDR_annotations <- Tg4510_FDR_annotations[row.names(new_common_sites_statsTg4510), ]
remove <- "NA"
new_Tg4510_FDR_annotations <- new_Tg4510_FDR_annotations[!(row.names(new_Tg4510_FDR_annotations) %in% remove), ]
write.csv(new_Tg4510_FDR_annotations, file = "/mnt/data1/isabel/RRBS/common_sites_statsTg4510.csv")

common_sites_statsJ20 <- J20_data[common_sites,]
new_common_sites_statsJ20 <- common_sites_statsJ20[ order(row.names(common_sites_statsJ20)), ]
new_J20_FDR_annotations <- J20_FDR_annotations[ order(row.names(J20_FDR_annotations)), ]
new_J20_FDR_annotations <- J20_FDR_annotations[row.names(new_common_sites_statsJ20), ]
remove <- "chr12:117177849"
new_J20_FDR_annotations <- new_J20_FDR_annotations[!(row.names(new_J20_FDR_annotations) %in% remove), ]
write.csv(new_J20_FDR_annotations, file = "/mnt/data1/isabel/RRBS/common_sites_statsJ20.csv")


# save genes for downstream analysis
# Annotated genes from common sites
genes_for_commonSites <- unique(new_Tg4510_FDR_annotations$Nearest_Gene)
# Common annotated genes
#common_genes
save(genes_for_commonSites, common_genes, file = "/mnt/data1/isabel/RRBS/CommonSitesandGenes.RData")