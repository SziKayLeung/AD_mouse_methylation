
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))
suppressMessages(library("stringr"))

DMRsGenotype <- read.csv("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMRs/rTg4510/DMRsGenotype.csv")
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_predictedMeth_rTg4510.RData") 
load(file = paste0(dirnames$figures, "PrnpPlots_corr_expression_methylation_byprobeGenotype.rData")); PrnpMethExp <- Prnpimage; rm(Prnpimage)


plotDMRGenotypeMedian <- function(predictedMeth = predictedMeth, DMR_df = DMR_df){
  
  library(BiSeq)
  betas <- as.matrix(methLevel(predictedMeth))
  coordinates <- as.data.frame(predictedMeth@rowRanges)
  pheno <- as.data.frame(predictedMeth@colData) %>% tibble::rownames_to_column(., var = "sample") %>% select(sample,Genotype,Age_months)
  DMR_df <- as.data.frame(DMR_df)
  coordinates$position <- rownames(coordinates)
  DMR_df$seqnames <- as.character(DMR_df$seqnames)
  
  
  plot_region <- function(region){
    
    DMR <- DMR_df[region,]
    
    # Select Cluster where this DMR lies
    #check cluster id
    Cluster_DMR_Chr <- coordinates[which(coordinates$seqnames == DMR[,"seqnames"]),]
    cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$start == DMR[,"start"]), "cluster.id"]
    
    #pull out cluster from coordinates and betas
    coordinates_cluster <- Cluster_DMR_Chr[which(Cluster_DMR_Chr$cluster.id == cluster),]
    betas_cluster <- betas[rownames(betas) %in% coordinates_cluster$position,]
    betas_coord_clus <- cbind(betas_cluster, coordinates_cluster)
    
    #add color to the pheno file for plotting
    
    if (!identical(pheno$sample, colnames(betas_cluster))) {
      stop("Rownames of betas and pheno matrices are not identical")
    }
    
    ## Plotting
    betas_coord_clus_melt <- betas_coord_clus %>% select(-seqnames, -end, -strand, -width, -cluster.id, -position) %>%
      reshape2::melt(., id = "start") 
    betas_coord_clus_melt <- na.omit(betas_coord_clus_melt)
    betas_coord_clus_melt <- betas_coord_clus_melt %>% mutate(methylation = as.numeric(value) * 100) %>% merge(., pheno, by.x = "variable", by.y = "sample")
    p <- ggplot(betas_coord_clus_melt, aes(x = start, y = methylation, colour = Genotype)) + geom_point() +
      theme_classic() + labs(x = unique(betas_coord_clus$seqnames), y = "Methylation (%)") +
      stat_summary(data=betas_coord_clus_melt, aes(x=start, y=methylation, group=Genotype), fun ="median", geom="line", linetype = "dotted", size = 1.5) +
      annotate("rect", xmin = DMR$start, xmax = DMR$end, ymin = 0, ymax = 100, fill = "gray", alpha = 0.2) +
      scale_colour_manual(values = c("black","#00AEC9"))
    
    return(p)
  }
  
  allplots <- lapply(DMR_df[["X"]],function(x) plot_region(x))
  
  return(allplots)
}


# Annotated DMP
annotatedDMP <- lapply(list.files(path = dirnames$annotated, full.names = T), function(x) read.csv(x))
names(annotatedDMP) <- lapply(list.files(path = dirnames$annotated), function(x) word(x,c(1),sep=fixed(".")))

# Ank1 reference gtf 
# subsetted Ank1 annotations from M22 reference GENCO
PrnprefGtf <- "/lustre/projects/Research_Project-MRC148213/sl693/reference/annotation/MousePrnp_gencodeM22.gtf"
PrnprefGtf <- as.data.frame(rtracklayer::import(PrnprefGtf))

# Table of Prnp 
PrnpDMP <- lapply(annotatedDMP, function(x) x %>% filter(ChIPseeker_GeneSymbol == "Prnp") %>% select(Position))
PrnpDMP <- bind_rows(PrnpDMP, .id = "file") %>% filter(grepl("rTg4510", file)) %>% mutate(seqnames = word(Position,c(1),sep=fixed(":")), 
                                                                                          start = as.numeric(word(Position,c(2),sep=fixed(":")))) %>% 
  mutate(end = start, transcript_id = word(file,c(3),sep=fixed("_")), strand = "+") %>% 
  select(transcript_id, start, end, seqnames, strand)



## ------------ p1: top ranked DMP -------

p1 <- create_df_plot("chr12:115283603") %>% 
  ggplot(., aes(x = Genotype, y = methylation, fill = Genotype)) + geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") +
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


## ------------ p2: Prnp DMR -------

plotDMR <- plotDMRGenotypeMedian(predictedMeth = predictedMeth, DMR_df = DMRsGenotype)
p2 <- plotDMR[[4]] + theme(legend.position = "None")


## ------------ p3: Prnp Pyro ------------

ank1PyroPos <- c(
  `Pos1Meth` = "chr2:131910162",
  `Pos2Meth` = "chr2:131910164",
  `Pos3Meth` = "chr2:131910180",
  `Pos4Meth` = "chr2:131910201"
)

p3 <- input_pyro$prnp %>% mutate(Sample.group = Group.ID) %>% select(Age, Sample.group, contains("Pos")) %>% 
  reshape2::melt(id = c("Age","Sample.group"), variable.name = "Position", value.name = "methylation") %>% 
  mutate(Sample.group = factor(Sample.group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = Sample.group, y = methylation, fill = Sample.group)) + geom_boxplot(aes(fill = Sample.group), outlier.shape = NA) + 
  facet_grid(~Position, labeller = as_labeller(ank1PyroPos)) +
  #geom_point(aes(fill = Sample.group), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") + theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


## ------------ p4: Prnp tracks ------------

# using ggtranscript for visaulilsation
# ENSMUST00000110688.8 representative ANK1 transcript (longest)
gexons <- PrnprefGtf  %>% dplyr::filter(type == "exon") %>% select(transcript_id, start, end, seqnames, strand) %>% 
  filter(transcript_id != "ENSMUST00000142070.1")
gexons <- rbind(gexons,PrnpDMP)


# plot tracks 
p4 <- ggtranscript_plot(gexons)

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

Prnplegend <-g_legend(PrnpMethExp[[1]])
PrnpMethExp <- lapply(PrnpMethExp, function(x) x + labs(x = NULL, y = NULL) + theme(legend.position = "None"))
PrnpMethExpListPlots <- plot_grid(plotlist = PrnpMethExp[1:6])
y.grob <- textGrob("RNA-Seq normalized counts", gp=gpar(fontsize=12), rot=90)
x.grob <- textGrob("Methylation (%)", gp=gpar(fontsize=12))
p5  <- grid.arrange(Prnplegend, arrangeGrob(PrnpMethExpListPlots, left = y.grob, bottom = x.grob), nrow = 2, heights=c(1, 10))


## ------------ Concatenated plot ------------

allP <- plot_grid(
  plot_grid(p1,p2,p3, labels = c("A","B","C"),scale = 0.95, nrow=1, rel_widths = c(0.2,0.4,0.4)),
  plot_grid(p4, labels = c("D"),scale = 0.95),
  plot_grid(p5,labels = c("E"), scale = 0.95),
  rel_heights = c(0.3,0.2,0.5),nrow=3)


pdf(paste0(dirnames$script,"Figures/MainFigures3.pdf"), width = 10, height = 15)
allP
dev.off()
