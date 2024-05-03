# title: "RRBS - BiSeq J20 - individual plots DMPs"
# author: "Isabel Castanho (I.S.Castanho@exeter.ac.uk)"

load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_predictedMeth_J20.RData") # predictedMeth

color_J20_TG <- "#FF5A62"


# Get phenotypic data
## Phenotypic data
coldata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_coldata_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)
levels(coldata$Age)

# ENTORHINAL CORTEX
class(coldata$ECX)

# Remove samples that have NA for pathology
coldata_ECX <- coldata[,c("Genotype", "Age_months", "Group_ID", "Histology_no", "ECX")]
coldata_ECX_clean <- na.omit(coldata_ECX)
coldata_pathology <- coldata_ECX_clean

# Prepare methylation the data
## get data
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_RRBSbetasComplete.RData")
# RRBS_completebetas
RRBS_completebetas_ECX <- RRBS_completebetas[,rownames(coldata_pathology)]

## Prepare data to plot
data <- RRBS_completebetas
coldata <- coldata
identical(rownames(coldata), colnames(data))


################################################# Genotype ################################################# 
## Remove sites with methylation difference < 0.02 (2%)
ResultsGenotype <- read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsInteractionModel_J20_sig_genotype.csv",
                            header = TRUE, row.names = 1)
colnames(ResultsGenotype)[1] <- "Location"
FilteredResultsGenotype <- ResultsGenotype[which(abs(ResultsGenotype[,"meth.diff.Genotype"]) > 0.02),]

## Get annotated tables
AnnotatedResultsGenotype <- read.csv("/lustre/projects/Research_Project-191406/EmmaW/Annotated/J20/DMPs/With_3000bpTSS_Annotations/DMPsInteractionModel_J20_sig_genotype_3000tssAnno.csv",
                                     header = TRUE)

AnnotatedResultsGenotype <- AnnotatedResultsGenotype[, -c(1:2)] # delete columns 1 through 2

FilteredAnnotatedResultsGenotype <- AnnotatedResultsGenotype[which(abs(AnnotatedResultsGenotype[,"meth.diff.Genotype"]) > 0.02),]

write.csv(FilteredAnnotatedResultsGenotype,
     file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/Annotated/J20/J20-DMPsGenotype_annotatedMethylDiffFiltered.csv")

## data to plot
DF <- FilteredAnnotatedResultsGenotype

## Top 500 FDR
# sort by FDR and then methylation difference
DF <- DF[
  with(DF, order(FDR_adj_genotype, -abs(meth.diff.Genotype))),
  ]
dataPlot <- DF[1:500,]
sites_to_plot <- dataPlot$Location
genes_to_plot <- as.character(dataPlot$biomart_Nearest_gene)

pdf(file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/IndividualPlots/J20/J20_GenotypeTop500FDR_MethylDiffFiltered.pdf", width=5, height=10)
  for (index in 1:length(sites_to_plot)) {
    site <- sites_to_plot[index]
    gene <- genes_to_plot[index]
    data_plot <- cbind(coldata[,c("Age_months","Genotype", "Group_ID")], t(data[as.character(site),]))
    data_plot <- data_plot[complete.cases(data_plot), ]
    #dev.new(width=5, height=10)
    par(mar=c(5,5,5,3))
    plot(y=data_plot[,as.character(site)], x=as.numeric(data_plot$Genotype),
         col = c("black", color_J20_TG)[data_plot$Genotype],
         pch = 19,
         cex=2,
         xlim=c(0.5,2.5),
         ylim=c(0, 1),
         axes=FALSE, ann=FALSE)
    axis(side=1, at=1:2, labels=c("WT", "TG"), cex.axis=2)
    axis(2, cex.axis=2)
    title(main = paste(gene, " ", site), ylab = "Methylation (%)", cex.main=2, cex.lab=2)
    boxplot(data_plot[,as.character(site)] ~ data_plot$Genotype,
    	axes=FALSE, ann=FALSE,
    	col = c("grey", color_J20_TG)[data_plot$Genotype],
    	add = TRUE, at = 1:2)
  }
dev.off()


################################################# Interaction (Genotype *Age) ################################################# 
## Remove sites with methylation difference < 0.02 (2%) between all WT timepoint 4 TG (all WT and 12 months TG)
ResultsInteraction <- read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsInteractionModel_J20_sig_interaction.csv",
                               header = TRUE, row.names = 1)
colnames(ResultsInteraction)[1] <- "Location"

## Get methylation mean difference data
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_DMPsInteraction_meanDiffTG-WT.RData")
# DatFrame
identical(rownames(DatFrame), rownames(ResultsInteraction)) # Must be TRUE

FilteredResultsInteraction <- ResultsInteraction[which(abs(DatFrame[,"meanDiffT4"]) > 0.02),]

## Get annotated tables
AnnotatedResultsInteraction <- read.csv("/lustre/projects/Research_Project-191406/EmmaW/Annotated/J20/DMPs/With_3000bpTSS_Annotations/DMPsInteractionModel_J20_sig_interaction_3000tssAnno.csv",
                                     header = TRUE)

AnnotatedResultsInteraction <- AnnotatedResultsInteraction[, -c(1:2)] # delete columns 1 through 2

identical(as.character(AnnotatedResultsInteraction$Location), rownames(DatFrame)) # Must be TRUE

AnnotatedResultsInteraction <- cbind(AnnotatedResultsInteraction, DatFrame$meanDiffT4)

colnames(AnnotatedResultsInteraction)[46] <- "meanDiffT4"

FilteredAnnotatedResultsInteraction <- AnnotatedResultsInteraction[which(abs(AnnotatedResultsInteraction[,"meanDiffT4"]) > 0.02),]

write.csv(FilteredAnnotatedResultsInteraction,
          file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/Annotated/J20/J20-DMPsInteraction_annotatedMethylDiffFiltered.csv")

## data to plot
DF <- FilteredAnnotatedResultsInteraction

## Top 500 FDR
# sort by FDR and then methylation difference
DF <- DF[
  with(DF, order(FDR_adj_interaction, -abs(meanDiffT4))),
  ]
dataPlot <- DF[1:500,]
sites_to_plot <- dataPlot$Location
genes_to_plot <- as.character(dataPlot$biomart_Nearest_gene)

pdf(file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/IndividualPlots/J20/J20_InteractionTop500FDR_MethylDiffFiltered.pdf", width=5, height=10)
for(index in 1:length(sites_to_plot)) {
  site <- sites_to_plot[index]
  gene <- genes_to_plot[index]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Group_ID")], t(data[as.character(site),]))
  data_plot <- data_plot[complete.cases(data_plot), ]
  means_WT <- c()
  means_TG <- c()
  for (age in c("6","8","10","12")) {
    means_WT <- c(means_WT, mean(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="WT"))[,as.character(site)]))
    means_TG <- c(means_TG, mean(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="TG"))[,as.character(site)]))
  }
  #dev.new(width=5, height=10)
  par(mar=c(5,5,5,3))
  plot(y=data_plot[,as.character(site)], x=as.numeric(data_plot$Age_months),
       col = c("black", color_J20_TG)[data_plot$Genotype],
       pch = 19,
       cex=2,
       xlim=c(0.5,4.5),
       ylim=c(0,1),
       axes=FALSE, ann=FALSE)
  axis(side=1, at=1:4, labels=c("6", "8", "10", "12"), cex.axis=2)
  axis(2, cex.axis=2)
  title(paste(gene, " ", site), xlab = "Age (months)", ylab = "DNA methylation", cex.main=2, cex.lab=2)
  lines(means_WT, lty=2, col="black", lwd=3)
  lines(means_TG, lty=2, col=color_J20_TG, lwd=3)
}
dev.off()


################################################# Pathology ################################################# 
## Prepare data to plot (keep only samples with pathology data)
data <- RRBS_completebetas_ECX
coldata <- coldata_pathology
identical(rownames(coldata_pathology), colnames(data))

## Remove sites with methylation difference < 0.02 (2%) between all WT timepoint 4 TG (all WT and 12 months TG)
ResultsPathology <- read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsPathology_J20_sig_pathology.csv",
                               header = TRUE, row.names = 1)
colnames(ResultsPathology)[1] <- "Location"

## Get methylation mean difference data
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_DMPsPathology_meanDiffTG-WT.RData")
# DatFrame

DatFrame <- DatFrame[rownames(ResultsPathology),]
identical(rownames(DatFrame), rownames(ResultsPathology)) # Must be TRUE

FilteredResultsPathology <- ResultsPathology[which(abs(DatFrame[,"meanDiffT4"]) > 0.02),]

## Get annotated tables
AnnotatedResultsPathology <- read.csv("/lustre/projects/Research_Project-191406/EmmaW/Annotated/J20/DMPs/With_3000bpTSS_Annotations/DMPsPathology_J20_sig_pathology_3000tssAnno.csv",
                                        header = TRUE)

AnnotatedResultsPathology <- AnnotatedResultsPathology[, -c(1:2)] # delete columns 1 through 2

identical(as.character(AnnotatedResultsPathology$Location), rownames(DatFrame)) # Must be TRUE

AnnotatedResultsPathology <- cbind(AnnotatedResultsPathology, DatFrame$meanDiffT4)

colnames(AnnotatedResultsPathology)[33] <- "meanDiffT4"

FilteredAnnotatedResultsPathology <- AnnotatedResultsPathology[which(abs(AnnotatedResultsPathology[,"meanDiffT4"]) > 0.02),]

write.csv(FilteredAnnotatedResultsPathology,
          file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/Annotated/J20/J20-DMPsPathology_annotatedMethylDiffFiltered.csv")

## data to plot
DF <- FilteredAnnotatedResultsPathology

## Top 500 FDR
# sort by FDR and then methylation difference
DF <- DF[
  with(DF, order(FDR_adj_pathology, -abs(meanDiffT4))),
  ]
dataPlot <- DF[1:500,]
sites_to_plot <- dataPlot$Location
genes_to_plot <- as.character(dataPlot$biomart_Nearest_gene)

pdf(file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/IndividualPlots/J20/J20_PathologyTop500FDR_MethylDiffFiltered.pdf", width=5, height=10)
for(index in 1:length(sites_to_plot)) {
  site <- sites_to_plot[index]
  gene <- genes_to_plot[index]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Group_ID")], t(data[as.character(site),]))
  data_plot <- data_plot[complete.cases(data_plot), ]
  means_WT <- c()
  means_TG <- c()
  for (age in c("6","8","10","12")) {
    means_WT <- c(means_WT, mean(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="WT"))[,as.character(site)]))
    means_TG <- c(means_TG, mean(subset(data_plot, subset=(data_plot$Age_months==age & data_plot$Genotype=="TG"))[,as.character(site)]))
  }
  #dev.new(width=5, height=10)
  par(mar=c(5,5,5,3))
  plot(y=data_plot[,as.character(site)], x=as.numeric(data_plot$Age_months),
       col = c("black", color_J20_TG)[data_plot$Genotype],
       pch = 19,
       cex=2,
       xlim=c(0.5,4.5),
       ylim=c(0,1),
       axes=FALSE, ann=FALSE)
  axis(side=1, at=1:4, labels=c("6", "8", "10", "12"), cex.axis=2)
  axis(2, cex.axis=2)
  title(paste(gene, " ", site), xlab = "Age (months)", ylab = "DNA methylation", cex.main=2, cex.lab=2)
  lines(means_WT, lty=2, col="black", lwd=3)
  lines(means_TG, lty=2, col=color_J20_TG, lwd=3)
}
dev.off()

### SCATTER PLOTS PATHOLOGY
pdf(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/ScatterPlots/J20/J20_PathologyTop500FDR_scatterPlots_Pathology.pdf", width=10, height=10)
  for (index in 1:length(sites_to_plot)) {
    site <- sites_to_plot[index]
    gene <- genes_to_plot[index]
    data_plot <- cbind(coldata[,c("Age_months","Genotype", "ECX")], t(data[as.character(site),]))
    data_plot <- data_plot[complete.cases(data_plot), ]
    #dev.new(width=10, height=10)
    par(mar=c(6,6,5,3), las=1,bty="l")
    plot(y=data_plot[,as.character(site)],
         x=data_plot$ECX,
         col = c("black", color_J20_TG)[data_plot$Genotype],
         #bg=c("white", color_J20_TG),
         pch = 19,
         cex=2,
         axes=FALSE,
         xlab = "Amyloid pathology in entorhinal cortex",
         ylab = "DNA methylation",
         xlim=c(0,2),
         ylim=c(0,1),
         mgp=c(4,1,0), cex.main=2, cex.lab=2)
    axis(side=1, cex.axis=2)
    axis(2, cex.axis=2)
    abline(fit <- lm(data_plot[,as.character(site)] ~ data_plot$ECX, data=data_plot), col='red')
    title(main = paste(gene, " ", as.character(site)), cex.main=2)
    mtext(coefs <- paste("R2=", signif(summary(fit)$r.squared, 2), ";", "P=", signif(summary(fit)$coefficients[1,4], 2)))
  }
dev.off()