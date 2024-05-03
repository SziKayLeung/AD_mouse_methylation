## Aisha Dahir A.N.Dahir@exeter.ac.uk
## Annotations to Horvath's coordinates are difficult with using biomart so will use chip package to help annotate

setwd("/mnt/data1/Dementia_Mouse/Array")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(tidyverse)
probe_mapping_file <- read_csv("/mnt/data1/Dementia_Mouse/MammalianArrayNormalizationTools/mappabilityData/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.1.csv", 
                               col_types = cols(.default = "c"))

#filter for mouse samples
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# peak can be granges object..
manifest <- species_probes
manifest$Chr <- gsub(":.*", "", manifest$MusMusculus)
manifest$Bp <- gsub(".*:", "", manifest$MusMusculus)
manifest <- manifest[-which(manifest$Chr %in% c("CHR_MG51_PATCH" ,
                                                "CHR_MG4200_PATCH", "CHR_MG3699_PATCH")),]

manifest$start <- manifest$Bp
manifest$end <- manifest$Bp
manifest$Chr <- paste("chr", manifest$Chr, sep = "")

rownames(manifest) <- manifest$probeID
manifest <- manifest[,-c(1:2)]
mm10_granges <- makeGRangesFromDataFrame(manifest,
                                         keep.extra.columns=T,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames","seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name"),
                                         start.field="start",
                                         end.field="end") #make this object to carry the liftover


peakAnno <- annotatePeak(mm10_granges,
                         TxDb=txdb, annoDb="TxDb.Mmusculus.UCSC.mm10.knownGene")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)

# Pull out the annotation 
peakAnno_df <- as.data.frame(peakAnno@anno)


# Convert the ensemble id startto ensemble name
EnsembleID <- peakAnno_df$transcriptId

library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(attributes = c('ensembl_transcript_id_version',
                               'ensembl_gene_id',
                               'external_transcript_name',
                               'external_gene_name'),
                filters = 'ensembl_transcript_id_version',
                values = EnsembleID,
                mart = mart)

write.csv(G_list, "Peak_transcript_Gene.csv")

peakAnno_df$Gene_Symbol <- rep(NA, nrow(peakAnno_df))

peakAnno_df$Gene_Symbol <- G_list$external_gene_name[match(peakAnno_df$transcriptId, G_list$ensembl_transcript_id_version)]


peakAnno_df$cpg <- rownames(peakAnno_df)
peakAnno_df <- peakAnno_df[,c(17,1:16)]

write.csv(peakAnno_df, "mm10_Manifest.csv")

# perhaps two genes to one sites??

peakAnno_df$Gene_Symbol2 <- rep(NA, nrow(peakAnno_df))

for(i in 1:nrow(peakAnno_df)){
  
  transcript <- peakAnno_df[i,"transcriptId"]
  
  genename <- G_list$external_gene_name[match(transcript, G_list$ensembl_transcript_id_version)]
  genename2 <- paste(unique(genename), collapse = " ")
  
  peakAnno_df[i,"Gene_Symbol2"] <- genename2
}


# The two gene name columns are not identical...
identical(peakAnno_df$Gene_Symbol, peakAnno_df$Gene_Symbol2)
peakAnno_df2 <- peakAnno_df[-which(peakAnno_df$Gene_Symbol %in% peakAnno_df$Gene_Symbol2),]