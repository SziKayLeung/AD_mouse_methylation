
##--- Annotations -------

suppressMessages(library(naniar))
suppressMessages(library(dplyr))

##-------------- ChipAnnotatePeaks

# Aim: original function called for chipseeker annotation used by E.Walker
# Input:
  # File = df for annotation (requiring position column = <chrA>:<Bposition>)  
# Output: 
  # peakAnno.df = df with column of genenames from chipseeker annotation

ChipAnnotatePeaks <- function(File,merge=NULL){
  
  #remove na
  #create GRanges object
  File <- File[!is.na(File$position),]
  seg <- GRanges(File$position)
  
  #add in feature info (promoter etc.) from chipseeker package
  peakAnno <- annotatePeak(seg,tssRegion = c(-1500, 1500),TxDb = txdb,
                           level = "transcript",assignGenomicAnnotation = TRUE,
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                           annoDb = "org.Mm.eg.db",
                           addFlankGeneInfo = FALSE, flankDistance = 5000, sameStrand = FALSE,ignoreOverlap = FALSE,ignoreUpstream = FALSE,
                           ignoreDownstream = FALSE,overlap = "TSS",verbose = TRUE)
  
  
  #extract more detailed feature info and add to output table
  peakAnno.df <- as.data.frame(peakAnno)
  peakAnno.df <- cbind(peakAnno.df, peakAnno@detailGenomicAnnotation)
  
  # create a new col for T/F only specified if genic OR within 1500bp window of tss
  peakAnno.df <- peakAnno.df %>% 
    mutate(InGeneBodyOr1500bpTSS = if_else(genic == T|(distanceToTSS >= -1500 & distanceToTSS <= 1500), T, F),
           position = paste0(seqnames,":", start))
  
  # create a new col for gene symbol only specified if genic OR within 1500bp window of tss
  peakAnno.df <- peakAnno.df %>% 
    mutate(InGeneBodyOr1500bpTSS_SYMBOL = if_else(InGeneBodyOr1500bpTSS == T, SYMBOL, ""))
  peakAnno.df <- peakAnno.df %>% replace_with_na(replace = list(InGeneBodyOr1500bpTSS_SYMBOL = ""))
  
  if(!is.null(merge)){
    message("Merging with the original df")
    commonCol <- intersect(colnames(File),colnames(peakAnno.df))
    File <- File %>% dplyr::select(-commonCol)
    peakAnno.df <- merge(File,peakAnno.df, by.x = merge, by.y = "position")
  }
  
  return(peakAnno.df)
}


##-------------- ChipMergeRefGtf

# Aim: Merge chipseeker annotation column with reference gtf using Ensembl transcript id 
# Input:
  # peakAnno.df = df: output from chipseeker annotation (ChipAnnotatePeaks)
  # refgtf = df: read in reference gtf
# Output: 
  # peakAnno.df_intronexon = subsetted annotations of peakAnno.df with additional column ("gene_id","gene_name") from merging with gtf
  # annotations are limited to probes in intron and exon

ChipMergeRefGtf <- function(peakAnno.df, refgtf){
  
  # extract "intron" and "exon" from annotation column
  peakAnno.df_intron <- peakAnno.df[grepl("intron",peakAnno.df$annotation),] 
  peakAnno.df_exon <- peakAnno.df[grepl("exon",peakAnno.df$annotation),] 
  
  # rbind both features and extract the ENSEMBL id from annotation 
  # E.g. Exon (ENSMUST00000034453.5/11459, exon 3 of 7) --> ENSMUST00000034453.5
  peakAnno.df_intronexon <- rbind(peakAnno.df_exon, peakAnno.df_intron) %>% 
    mutate(transcriptID_ref  = word(word(word(annotation,c(2), sep = fixed(" ")),c(1),sep = fixed("/")),c(2), sep = fixed("(")))
  
  # extract the transcript_id, gene_id, and gene_name from reference gtf
  genes = unique(refgtf[ ,c("transcript_id","gene_id","gene_name")])
  
  # merge 
  peakAnno.df_intronexon <- merge(peakAnno.df_intronexon, genes, by.x = "transcriptID_ref", by.y = "transcript_id", all.x = T)
  
  return(peakAnno.df_intronexon)
}


##-------------- ChipFinaliseAnno

# Aim: Merge output from ChipAnnotatePeaks with output from ChipMergeRefGtf
# i.e remove redundant probes in ChipAnnotatePeaks output df that are further annotated using ChipMergeRefGtf
# Input:
  # peakAnno.df = df: output from chipseeker annotation (ChipAnnotatePeaks)
  # peakAnno.df_intronexon = output from of ChipMergeRefGtf
  # stat.df = df: reported statistics with "position" column
# Output: 
  # mergedPeakAnno.df = finalised annotation table with chipseeker columns and additional annotation (gene_id, gene_name) from merging with reference
# NB: 
  # warning : 'select()' returned 1:many mapping between keys and columns, from multiple entries in the org.Mm.eg.db annotation that have the same Entrez-gene id --> not to worry as chipseeker picks the first one (https://support.bioconductor.org/p/100899/)
  
ChipFinaliseAnno <- function(peakAnno.df, peakAnno.df_intronexon, stat.df){
  
  cols2keep <- c("position","annotation","gene_id","gene_name","transcriptID_ref", "ENSEMBL","SYMBOL","GENENAME","transcriptId","distanceToTSS",
                 "genic","Intergenic","Promoter","fiveUTR",
                 "threeUTR","Exon","Intron","downstream","distal_intergenic","InGeneBodyOr1500bpTSS","InGeneBodyOr1500bpTSS_SYMBOL")
  
  # create <gene_id> and <gene_name> column for later rbind, given these two columns are in peakAnno.df_intronexon from merging with gtf
  # remove duplicated entries of probes that are further annotated in intron and exon from merging with gtf
  # note peakAnno.df_intronexon is the same df as peakAnno.df for those intron&exon probes, aside from two additional column (gene_id, gene_name)
  peakAnno.df <- peakAnno.df  %>% mutate(gene_id = "NA",gene_name="NA", transcriptID_ref="NA") %>% filter(!position %in% peakAnno.df_intronexon$position)
  
  # rbind the two dataframes, keeping only the columns of interest
  mergedPeakAnno.df <- rbind(peakAnno.df[cols2keep], peakAnno.df_intronexon[cols2keep])
  
  # merge annotated table with statistics table using "position" column
  mergedPeakAnno.df <- merge(stat.df, mergedPeakAnno.df, by = "position", all = T)
  
  # rename column names
  mergedPeakAnno.df <- mergedPeakAnno.df %>% 
    rename(
      Position = position,
      ChIPseeker_Annotation = annotation,
      RefGtf_GeneEnsembl = gene_id,
      RefGtf_GeneSymbol = gene_name, 
      RefGtf_TransEnsembl = transcriptID_ref,
      ChIPseeker_GeneEnsembl = ENSEMBL,
      ChIPseeker_GeneSymbol = SYMBOL,
      ChIPseeker_GeneName = GENENAME,
      ChIPseeker_TransEnsembl = transcriptId
    )
  
  # create column to if gene symbols are matching across both methods
  mergedPeakAnno.df <- mergedPeakAnno.df %>% mutate(
    Consensus_GeneSymbol = ifelse(as.character(RefGtf_GeneSymbol) == as.character(ChIPseeker_GeneSymbol),"TRUE","FALSE"),
    Consensus_GeneSymbol = ifelse(RefGtf_GeneEnsembl == "NA","NA",Consensus_GeneSymbol))
  

  
  return(mergedPeakAnno.df)
}


##-------------- ChipClassifyAnno

# Aim: Use chipseeker annotation column output to output simplified category
# category: 3'UTR, Distal Intergenic, Exon, Intron, Promoter, Downstream, Distal Intergenic
# output entry same as input intry except in string input containing <Intron> <Exon>
# Input:
  # annotation = str <3'UTR, Distal Intergenic, Intron (ENSMUSTXXX, intron X of X), Exon (ENSMUSTXXX, intron X of X), Downstream 2-3kb, ...>
# Output:
  # annotation = str <3'UTR, Distal Intergenic, Intron, Exon, Downstream, ...>

ChipClassifyAnno <- function(annotation){
  if(annotation %in% c("3' UTR","5' UTR","Distal Intergenic")){ 
    return(annotation) 
  } else if(grepl("Intron",annotation)) { 
    #pos = word(word(annotation, c(2), sep = fixed(",")),c(3),sep = fixed(" "))
    #ifelse(pos == "1", return("1st_intron"), return("other_intron")) 
    return("Intron")
  } else if (grepl("Exon",annotation)){
    #pos = word(word(annotation, c(2), sep = fixed(",")),c(3),sep = fixed(" "))
    #ifelse(pos == "1", return("1st_exon"), return("other_exon")) 
    return("Exon")
  } else if (annotation == "Promoter"){ 
    return("Promoter") 
  } else if (grepl("Downstream",annotation)){
    #pos = word(word(annotation, c(2), sep = fixed(",")),c(3),sep = fixed(" "))
    #ifelse(pos == "1", return("1st_exon"), return("other_exon")) 
    return("Downstream (<3kb)")
  } else if (annotation == "Promoter"){ 
    return("Promoter") 
  } else { 
    return(annotation) 
  }
  
}


##-------------- ChipClassifyAnno_df

# Aim: apply ChipClassifyAnno() to a chipseeker output for downstream plotting (ChipPlotAnno())
# Input:
  # df = output from ChipAnnotatePeaks() after running chipseeker; contains annotation column 
# Output:
  # subsetted input df <position, annotation, grpannotation> 
    # position = probe location 
    # annotation = original annotation from chipseeker output 
    # grpannotation = output from running ChipClassifyAnno()

ChipClassifyAnno_df <- function(df){
  
  df$grpannotation <- lapply(df$annotation, function(x) ChipClassifyAnno (x))
  df$grpannotation <- as.factor(as.character(df$grpannotation ))
  df <- df %>% dplyr::select(position, annotation, grpannotation )
  
  return(df)
}


##-------------- ChipPlotAnno

# Aim: plot chipseeker annotations of all and significant sites in rTg4510 and J20 mouse model
# Input:
  # ChipPlotAnnoOut = list of rTg4510_all, J20_all, J20_rrbs_genotype, rTg4510_rrbs_genotype
  # in each list ==> position; annotation; grpannotation  (output from running ChipPlotAnno_df)
# Output:
  # bar-plot of the percentage annotations of all RRBS sites and significant sites in rTg4510 and J20 mouse model

ChipPlotAnno <- function(ChipPlotAnnoOut, model=NULL){
  
  if(is.null(model)){
    # merge into one dataframe
    all <- rbind(ChipPlotAnnoOut$rTg4510_all %>% mutate(model = "rTg4510", sites = "all"),
                 ChipPlotAnnoOut$J20_all %>% mutate(model = "J20", sites = "all"),
                 ChipPlotAnnoOut$J20_rrbs_genotype  %>% mutate(model = "J20", sites = "DMP"),
                 ChipPlotAnnoOut$rTg4510_rrbs_genotype  %>% mutate(model = "rTg4510", sites = "DMP"))
  }else{
    all <- rbind(ChipPlotAnnoOut$rTg4510_all %>% mutate(model = "rTg4510", sites = "all"),
                 ChipPlotAnnoOut$rTg4510_rrbs_genotype  %>% mutate(model = "rTg4510", sites = "DMP"))
  }

  
  # datawrangle for plotting
  # tally the number of sites by model and annotation
  # calculate the percentage for each group
  # set levels for plot
  all_tally <- all %>% group_by(model, sites, grpannotation ) %>% tally() %>% 
    mutate(perc = n/sum(n) * 100) %>%
    mutate(model = factor(model, levels = c("rTg4510", "J20")),
           sites = factor(sites, levels = c("all", "DMP"), labels = c("All sites", "Significant sites")),
           grpannotation  = factor(grpannotation , 
                                     levels = c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR","Downstream (<3kb)", "Distal Intergenic")))
  
  # plot
  col_palette = c(wes_palette("Darjeeling2")[4], wes_palette("Rushmore1")[4], wes_palette("Darjeeling1")[2],
                  wes_palette("Royal1")[1], wes_palette("Darjeeling2")[1], wes_palette("Royal1")[4],wes_palette("Royal1")[2])
  if(is.null(model)){
    p <- ggplot(all_tally, aes(x = model, y = perc, fill = forcats::fct_rev(grpannotation ))) +  
      geom_bar(position = "stack", stat = "identity") + facet_grid(~sites) +
      labs(x = "Mouse  model", y = "Percentage of sites") + 
      scale_fill_manual(name = "Annotation", values = col_palette, guide = guide_legend(reverse = TRUE)) + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "grey50"),
            strip.background =element_rect(fill="white"))
  }else{
    p <- ggplot(all_tally, aes(x = sites, y = perc, fill = forcats::fct_rev(grpannotation ))) +  
      geom_bar(position = "stack", stat = "identity") +
      labs(x = NULL,y = "Percentage of sites") + 
      scale_fill_manual(name = "Annotation", values = col_palette, guide = guide_legend(reverse = TRUE)) +
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "grey50"),
            strip.background =element_rect(fill="white"))
  }
  
  return(p)
}
