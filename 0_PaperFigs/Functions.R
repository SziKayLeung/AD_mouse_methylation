suppressMessages(library(qqman))
suppressMessages(library(cowplot))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggrepel))
suppressMessages(library(ggh4x))

pastelColours <- brewer.pal(4, "Pastel2")


color_Tg4510_TG <- "#00AEC9"

# is one list in another list
VectorIntersect <- function(v,z) {
  unlist(lapply(unique(v[v%in%z]), function(x) rep(x,min(sum(v==x),sum(z==x)))))
}
is.contained <- function(v,z) {length(VectorIntersect(v,z))==length(v)}

# heatmap
suppressMessages(library(pheatmap))

# gene tracks
suppressMessages(library(ggbio))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

label_colour <- function(var){
  if(var %in% c("Tg4510","rTg4510")){colour = "#00AEC9"}else{
    if(var == "J20"){colour = "#FF5A62"}else{
    }}
  return(colour)
}

# plot RRBS CpG sites by annotation
plot_annotate_sites <- function(){
  # simplify annotations from chipseeker to exon and intron
  anno_rrbs_all <- lapply(anno_rrbs_all, function(x) x %>% 
                            mutate(annotation_simple = ifelse(grepl("Exon", annotation), "Exon", annotation),
                                   annotation_simple = ifelse(grepl("Intron", annotation), "Exon", annotation_simple)))
  
  # tally and percentage
  anno_stats <- lapply(anno_rrbs_all, function(x) as.data.frame(x %>% group_by(annotation_simple) %>%
                                                                  summarise(count = n()) %>% 
                                                                  mutate(perc = count/sum(count) * 100)))
  
  # binomial test, number of promoter sites
  rTg4510_promoter = anno_stats$rTg4510[anno_stats$rTg4510$annotation_simple == "Promoter","count"]
  J20_promoter = anno_stats$J20[anno_stats$J20$annotation_simple == "Promoter","count"]
  binom.test(rTg4510_promoter, sum(anno_stats$rTg4510$count), p = 0.5, alternative = "two.sided")
  binom.test(J20_promoter, sum(anno_stats$J20$count), p = 0.5, alternative = "two.sided")
  
  p <- bind_rows(anno_stats$rTg4510 %>% mutate(model = "rTg4510"), anno_stats$J20 %>% mutate(model = "J20")) %>%
    ggplot(., aes(x = model, y = perc, fill = annotation_simple)) + geom_bar(stat = "identity") +
    labs(x = "Mouse model", y = "Percentage of CpG sites (%)") + theme_classic() +
    scale_fill_discrete(name = "Annotations")  
  
  return(p)
}


# hierarchal clustering by the top 1000 most differentially expressed probes
cluster_DMP <- function(model, bothBeta = NULL, rrbsBeta = NULL, arrayBeta = NULL, phenotypeInput, lstPositions, clusterNum=1000){
  
  colourPoints <- label_colour(model)
  
  if(!is.null(bothBeta)){
    m <- bothBeta %>% filter(row.names(.) %in% lstPositions[1:clusterNum])
    
  }else{
    common_samples <- intersect(colnames(rrbsBeta), colnames(arrayBeta))
    
    m <- rbind(arrayBeta %>% filter(row.names(.) %in% lstPositions[1:clusterNum]) %>% dplyr::select(all_of(common_samples)),
               rrbsBeta %>% filter(row.names(.) %in% lstPositions[1:clusterNum]) %>% dplyr::select(all_of(common_samples)))
  }
 
  
  ann_colors = list(
    Genotype = c(WT = "black", TG = colourPoints)
  )
  
  p <- pheatmap(m, annotation = phenotypeInput[,c("Genotype","Age_months")],
                show_rownames = FALSE, show_colnames = FALSE,
                annotation_colors = ann_colors)
  
  return(p)
}

## ------ manhattan plots -----

prepare_manhattan <- function(df,p_val,type){
  
  if(length(grep("location",names(df),value=TRUE))==1){
    df <- df %>% mutate(Location = location)
  }
  
  if(length(grep("Position",names(df),value=TRUE))==1){
    df <- df %>% mutate(Location = Position)
  }
  
  if(length(grep("position",names(df),value=TRUE))==1){
    df <- df %>% mutate(Location = position)
  }
  
  result <- df %>% 
    mutate(SNP = Location, 
           CHR = stringr::str_remove(word(Location,c(1),sep = stringr::fixed(":")),"chr"),
           BP = stringr::word(Location,c(2),sep = stringr::fixed(":"))) %>% 
    # patch ones ("CHR_MG51_PATCH", "CHR_MG4200_PATCH","CHR_MG3699_PATCH")
    filter(!is.na(BP)) %>%
    # Remove chrY (all samples are females) and chrM
    filter(!CHR %in% c("X","Y","M")) %>% 
    mutate(CHR = as.numeric(CHR),BP = as.numeric(BP))
  
  if (type == "FDR_correct"){
    result$mFDR  <- p.adjust(result[,p_val], method = "fdr")
  } else if (type == "FDR_present"){
    result$mFDR  <- result[[p_val]]
  }else{
    NULL
  }
  
  return(result)
}

plot_manhattan_final <- function(res, term){
  
  pvalueTerm <- ifelse(term == "Genotype", "P.value_Genotype", "P.value_Pathology")
  prep_mhat = prepare_manhattan(res[[term]], pvalueTerm, "FDR_present")
  cumu <- cumulative(prep_mhat)
  
  # parameters
  
  # label gene names
  cumu <- cumu %>% mutate(anno = ifelse(cumu[[pvalueTerm]] < gwas_sig, as.character(ChIPseeker_GeneSymbol),""),
                          alpha = as.factor(ifelse(cumu[[pvalueTerm]] < gwas_sig, TRUE, FALSE)))
  
  # basic plot
  x.var <- rlang::sym(quo_name(enquo(pvalueTerm)))
  plot <- ggplot(cumu, aes(x = bp_cum, y = -log10(!!x.var), color = Platform, label = anno)) +
    geom_point(aes(shape = Platform, alpha = alpha, size = 2)) +
    scale_shape_manual(values = c(16, 17)) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) 
  
  # axis formats
  axis_set <- cumu %>% group_by(CHR) %>% summarize(center = mean(bp_cum))
  axis_edge <- data.frame(cumu %>% group_by(CHR) %>% summarise(edge = max(bp_cum)))
  ylim <- abs(floor(log10(min(cumu[[pvalueTerm]], na.rm = T)))) + 2
  
  # finalise plot
  plot <- plot +  
    geom_label_repel(aes(label=anno), size=5, show_guide  = FALSE, box.padding = 0.5, max.overlaps = Inf) + 
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0)) +
    #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
    geom_hline(yintercept = -log10(gwas_sig), color = "black", linetype = "dashed") +
    geom_vline(xintercept = c(axis_edge$edge),linetype = "dotted", colour = "lightgrey") +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = expression(-log[10](italic(p)))) + 
    theme_classic() +
    theme( 
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 8, vjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) + guides(alpha="none",text="none",size="none") 
  
  return(plot)
  
}

datawrangle_manhattan <- function(dat, pvalueCol, platform){
  dat <- dat %>% dplyr::select(Position,  pvalueCol) 
  dat$CHR <- stringr::word(dat$Position,c(1),sep=stringr::fixed(":"))
  dat$BP <- stringr::word(dat$Position,c(2),sep=stringr::fixed(":"))
  dat[,"CHR"][which(dat[,"CHR"] == "chrX")]<-"chr23"
  dat[,"CHR"][which(dat[,"CHR"] == "chrY")]<-"chr24"
  dat$CHR <- as.numeric(str_remove(dat$CHR,"chr"))
  dat$BP <- as.numeric(dat$BP)
  dat <- dat %>% mutate(Platform = platform)
  return(dat)
}

plot_merged_manhattan <- function(mergedManhattanPlot, annoRRBS, annoArray, term){
  
  # remove NAs
  mergedManhattanPlot <- mergedManhattanPlot[!is.na(mergedManhattanPlot$CHR),]
  
  # cumulative position of SNPs
  don <- mergedManhattanPlot %>% 
    
    # Compute chromosome size
    dplyr::group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(mergedManhattanPlot, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  
  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # gene annotations
  donRRBS <- merge(don[don$Platform == "RRBS",], 
                   annoRRBS[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  donArray <- merge(don[don$Platform == "Array",], 
                    annoArray[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  don <- rbind(donRRBS, donArray)
  
  if(term != "Genotype"){
    don <- don %>% dplyr::arrange(p.val.Pathology) 
  }else{
    don <- don %>% dplyr::arrange(p.val.Genotype)
  }
  top <- unique(don$ChIPseeker_GeneSymbol)[1:200]
  don <- don %>% mutate(label = ifelse(ChIPseeker_GeneSymbol %in% top, ChIPseeker_GeneSymbol, NA))
  don <- don %>% mutate(label_unique = if_else(duplicated(label), NA_character_, label)) 
  
  # plot
  if(term != "Genotype"){
    p <- ggplot(don, aes(x=BPcum, y = -log10(p.val.Pathology)))
  }else{
    p <- ggplot(don, aes(x=BPcum, y = -log10(p.val.Genotype)))
  }
  
  p <- p +
    
    # Show all points
    geom_point( aes(color=as.factor(Platform)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("red", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    geom_text_repel(data = don[!is.na(don$label_unique), ], 
                    aes(label = label_unique), 
                    size = 3, 
                    box.padding = 0.3, 
                    point.padding = 0.5) +
    
    # Custom the theme:
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))))
  
  return(p)
}

plotGeneTrackDMP <- function(sigResults, betaMatrix, phenotypeFile, gene, transcript, boxplot = FALSE, colour = FALSE,
                             pathology = FALSE, position = NULL){
  
  if(isFALSE(colour)){
    colourbox = "yellow"
  }else{
    colourbox <- label_colour(colour)
  }
  
  # extract positions from beta matrix
  if(is.null(position)){
    dat <- betaMatrix %>% filter(row.names(betaMatrix) %in% sigResults[sigResults$ChIPseeker_GeneSymbol %in% gene,"Position"])
  }else{
    dat <- betaMatrix %>% filter(row.names(betaMatrix) %in% position)
  }

  # split to get the coordinates from the position <chrX:YY>
  dat <- dat %>% tibble::rownames_to_column(., var = "position") %>% reshape2::melt(variable.name = "sample",value.name = "methylation", id = "position")
  dat <- merge(dat, phenotypeFile, by.y = 0, by.x = "sample")
  dat$coordinate <- stringr::str_split_i(dat$position,":",2)
  dat$chr <- stringr::str_split_i(dat$position,":",1)
  
  
  # extract the transcript of interest from txdb
  gr <- subset(transcripts(txdb), tx_name == transcript)
  grdf <- as.data.frame(gr)
  
  
  # gene track (note reduce: collapsed all the exons within that vicinity from transcript)
  # stat = "reduce"
  gene_track <- ggplot() + 
    geom_alignment(TxDb.Mmusculus.UCSC.mm10.knownGene, which = gr, label = FALSE) + 
    theme_bw() + 
    labs(subtitle = gene) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank(),
          plot.subtitle = element_text(face = "italic")) 

  
  # min-value and max-value from the DMP range
  minvalue = min(dat$coordinate)
  maxvalue = max(dat$coordinate)
  
  # box the DMP region
  gene_track <- gene_track +
    geom_rect(data = as.data.frame(grdf), aes(xmin = as.numeric(minvalue) , xmax = as.numeric(maxvalue), ymin = -Inf, ymax = Inf), 
              fill = colourbox, alpha = 0.3, 
              colour = colourbox)
  
  if(isFALSE(boxplot)){
    
    p <- ggplot(dat, aes(x = as.numeric(coordinate), y = methylation, colour = Genotype)) +
      geom_point() +
      scale_color_manual(values=c("black", colourbox)) +
      theme_classic() +
      stat_summary(aes(col = Genotype, group = Genotype), fun.y = mean, geom = "smooth", linetype = "dotted", colour = colourbox) +
      theme_classic() + 
      theme(panel.border = element_rect(fill = NA, color = "white", linetype = "dotted"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank()) +
      labs(y = "Methylation", x = paste0("Co-ordinate (", dat$chr[1],")")) +
      theme(legend.position = "None", 
            #panel.background = element_rect(colour = colourbox, fill = alpha("white",0.1))
      ) 
    
    
  }else{
    if(isFALSE(pathology)){
      p <- plot_DMP(betaMatrix, phenotypeFile, position = unique(as.character(dat$position)), pathology = FALSE, model = colour) 
    }else{
      p <- plot_DMP(betaMatrix, phenotypeFile, position = unique(as.character(dat$position)), pathology = TRUE, model = colour) 
    }
  }
  
  output <- plot_grid(gene_track,p,nrow=2, rel_heights = c(0.3,0.7))
  return(output)
}


## ------ correlate pyrosequencing and rrbs -----

plot_pyro_rrbs_corr <- function(inputPyro, inputPyroPos, rrbsBeta, inputPhenotype){
  
  color_Tg4510_TG <- "#00AEC9"
  
  # data-wrangle input pyrosequencing data 
  dat <- inputPyro %>% mutate(Sample.group = Group.ID) %>% dplyr::select(SAMPLE, contains("Pos")) %>% 
    reshape2::melt(id = c("SAMPLE"), variable.name = "prnpPosition", value.name = "methylationPyro") %>% 
    merge(., inputPyroPos, by = "prnpPosition") %>% 
    dplyr::select(SAMPLE, methylationPyro, Position) %>%
    mutate(SamplePos = paste0(SAMPLE,Position))
  
  # data-wrangle rrbsBeta dataframe
  dat2 <- rrbsBeta %>% filter(row.names(.) %in% dat$Position) %>% tibble::rownames_to_column(., var = "Position") %>% 
    reshape2::melt(variable.name = "SAMPLE", value.name = "RRBSmethylation") %>%
    mutate(SamplePos = paste0(SAMPLE,Position))
  
  merged <- merge(dat, dat2, by = "SamplePos")
  merged <- merge(merged, inputPhenotype %>% tibble::rownames_to_column(., var = "SAMPLE.x"), by = "SAMPLE.x")
  
  # mutate to RRBS percentage
  merged <- merged %>% mutate(RRBSmethylation = RRBSmethylation * 100)
  #print(merged)
  
  output <- list()
  for(i in 1:length(unique(merged$Position.x))){
    pos <- unique(merged$Position.x)[i]
    dat <- merged[merged$Position.x == pos,]
    message("position: ", pos)
    print(cor.test(dat$methylationPyro, dat$RRBSmethylation))
  #  output[[i]] <- ggplot(dat, aes(x = RRBSmethylation, y = methylationPyro, colour = Genotype, group = Genotype)) + geom_point(size = 3) + 
  #    theme_classic() +
  #    labs(y = "Pyrosequencing methylation (%)", x = "RRBS methylation (%)", subtitle = pos) +
  #    theme(strip.background = element_blank()) +
  #    geom_smooth(method=lm, formula = y~poly(x,3),fill = "white", linetype = "dotted") +
  #    scale_color_manual(values=c("black", color_Tg4510_TG)) 
  }
  
  p <- ggplot(merged, aes(x = RRBSmethylation, y = methylationPyro, colour = Genotype, group = Genotype)) + geom_point(size = 3) + 
    theme_classic() +
    labs(y = "Pyrosequencing methylation (%)", x = "RRBS methylation (%)") +
    theme(strip.background = element_blank()) +
    geom_smooth(method=lm, formula = y~poly(x,3),fill = "white", linetype = "dotted") +
    scale_color_manual(values=c("black", color_Tg4510_TG)) +
    facet_grid(~Position.x)
 
  
  return(p)
}


## ------ epigenetic clock -----

plot_clock <- function(clock, tissue, model, boxplot = TRUE){
  

  y.var <- sym("DNAmAgeClockCortex")
  y.lab <- "DNAm Age Clock Cortex"

  
  if(model %in% c("rTg4510","Tg4510")){
    colour <- label_colour("Tg4510")
  }else{
    colour <- label_colour("J20")
  }
  
  if(!isFALSE(boxplot)){
    p <- ggplot(clock, aes(x = as.factor(Age_months), y = !! y.var, colour = Genotype)) + 
      geom_boxplot(outliers = FALSE) +
      geom_point(position=position_jitterdodge()) +
      labs(x = "Age (Months)", y = y.lab, colour = "Genotype") +
      theme_classic() +
      scale_colour_manual(values = c("black",colour)) 
  }else{
    p <- ggplot(clock, aes(x = Age_months, y = DNAmAgeClockCortex, colour = Genotype)) + geom_point() +
      theme_classic() + labs(x = "Age", y = "DNAm Age Clock Cortex", colour = "Genotype") +
      geom_smooth(method='lm', formula= y~x, se = FALSE) +
      scale_colour_manual(values = c("black",colour)) 
  }

  
  return(p)
  
}


clock_acceleration <- function(clock, mouse, tissue){
  
  if(mouse %in% c("rTg4510","Tg4510")){
    colour <- label_colour("Tg4510")
  }else{
    colour <- label_colour("J20")
  }
  
  if(tissue == "ECX"){
    x.var <- sym("Pathology_ECX")
    dat <- merge(clock, phenotype_path[[mouse]], by.y = "Sample_ID_ECX", by.x = "SampleID")
  }else{
    x.var <- sym("Pathology_HIP")
    dat <- merge(clock, phenotype_path[[mouse]], by.y = "Sample_ID_HIP", by.x = "SampleID")
  }
  
  p <- dat %>% mutate(rate =  DNAmAgeClockCortex/Age) %>% 
    ggplot(., aes(x = !! sym(x.var), y = rate, colour = Genotype.y)) + geom_point() +
    labs(x = paste0("Pathology in ", tissue), y = "Acceleration age", colour = "Genotype") +
    theme_classic() +
    scale_colour_manual(values = c(colour, "black")) 
  
  return(p)
}


clock_stats <- function(clock, age, modelTissue){
  print(paste0("Clock on ", modelTissue, " at ", age, " months"))
  print(with(clock %>% filter(Age_months == age), shapiro.test(DNAmAgeClockCortex[Genotype == "TG"])))
  print(with(clock %>% filter(Age_months == age), shapiro.test(DNAmAgeClockCortex[Genotype == "WT"])))
  print(t.test(DNAmAgeClockCortex ~ Genotype, data = clock %>% filter(Age_months == age)))
}



## ------ GO -----

extract_postions_as_bed <- function(Position, path){
  
  dat <- data.frame(Position)
  colnames(dat) <- "Position"
  dat <- dat %>% mutate(chr = word(Position,c(1),sep=fixed(":")), 
                        pos = word(Position,c(2),sep=fixed(":")),
                        pos2 = pos) %>% 
    dplyr::select(chr, pos, pos2)
  
  write.table(dat, path, col.names = F, row.names = F, quote = F, sep = "\t")
  
}
