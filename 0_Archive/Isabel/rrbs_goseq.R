#Isabel Castanho I.Castanho@exeter.ac.uk

# Functional annotation and gene ontology | goseq from common sites and genes from RRBS genotype analysis (t-test)

# R3.4.3
setwd("/mnt/data1/isabel/RRBS/")

## SETUP
source("http://bioconductor.org/biocLite.R") 

library(DESeq2) # DESeq2_1.16.1
library(GO.db)
library(goseq) # goseq_1.30.0 
library (org.Mm.eg.db)
library(annotate)
library(geneLenDataBase) # there is no package called ‘geneLenDataBase’

biocLite("Organism.dplyr")
library("Organism.dplyr")

load(file = "/mnt/data1/isabel/RRBS/CommonSitesandGenes.RData")
write.csv(genes_for_commonSites, file = "/mnt/data1/isabel/RRBS/genes_for_commonSites.csv")
write.csv(common_genes, file = "/mnt/data1/isabel/RRBS/commonGenes.csv")


### Function to run GOseq in overlapping genes from rTg4510 data and J20 data


# To see which genome/gene identifier combinations are in the local database
supportedOrganisms() ########################## ERROR ##########################

supportedOrganisms()[supportedOrganisms()$Genome=="mm10",]

# Fitting the Probability Weighting Function (PWF)
pwf=nullp(genes_vector,"mm10","geneSymbol")
head(pwf)




# Annotated genes for 364 common DMPs
#genes_for_commonSites

# Common annotaded genes (5007) 
#common_genes