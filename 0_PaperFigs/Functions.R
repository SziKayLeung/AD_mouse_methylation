suppressMessages(library(qqman))
suppressMessages(library(cowplot))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggrepel))
suppressMessages(library(ggh4x))
suppressMessages(library(extrafont))
suppressMessages(library(showtext))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap)) # heatmap
# gene tracks
suppressMessages(library(ggbio))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

pastelColours <- brewer.pal(4, "Pastel2")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=16),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 12),
                 axis.text.x= element_text(size=16),
                 axis.text.y= element_text(size=16),
                 plot.title = element_text(size=16),
                 plot.subtitle = element_text(size=16))


color_Tg4510_TG <- "#00AEC9"

# is one list in another list
VectorIntersect <- function(v,z) {
  unlist(lapply(unique(v[v%in%z]), function(x) rep(x,min(sum(v==x),sum(z==x)))))
}
is.contained <- function(v,z) {length(VectorIntersect(v,z))==length(v)}


label_colour <- function(var){
  if(var %in% c("Tg4510","rTg4510")){colour = "#00AEC9"}else{
    if(var == "J20"){colour = "#FF5A62"}else{
    }}
  return(colour)
}

## ------ plot sites -----

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


plot_gene_track <- function(betaMatrix, phenotypeFile, position, colour, gene, transcript){
  
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
          text = element_text(size = 16),
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
  
  return(gene_track)
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
          text = element_text(size = 16),
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
      stat_summary(aes(colour = Genotype, group = Genotype), fun.y = mean, geom = "smooth", linetype = "dotted") +
      theme_classic() + 
      theme(panel.border = element_rect(fill = NA, color = "white", linetype = "dotted"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank()) +
      labs(y = "Methylation", x = paste0("Co-ordinate (", dat$chr[1],")")) +
      mytheme +
      theme(legend.position = "None", 
            #panel.background = element_rect(colour = colourbox, fill = alpha("white",0.1))
      ) 
    
    
  }else{
    if(isFALSE(pathology)){
      p <- plot_DMP(betaMatrix, phenotypeFile, position = unique(as.character(dat$position)), pathology = FALSE, model = colour) + mytheme
    }else{
      p <- plot_DMP(betaMatrix, phenotypeFile, position = unique(as.character(dat$position)), pathology = TRUE, model = colour) + mytheme
    }
  }
  
  output <- plot_grid(gene_track,p,nrow=2, rel_heights = c(0.3,0.7))
  return(output)
}

plot_DMP_byTissue <- function(ECXbetaMatrix, HIPbetaMatrix, ECXphenotypeFile, HIPphenotypeFile, position, 
                              interaction = FALSE, pathology = FALSE, model = "rTg4510", gene = NULL, transcript = NULL){
  
  ECX_dat <- merge_beta_phenotype(ECXbetaMatrix, ECXphenotypeFile, position) %>% mutate(tissue = "ECX")
  HIP_dat <- merge_beta_phenotype(HIPbetaMatrix, HIPphenotypeFile, position)%>% mutate(tissue = "HIP")
  dat <- rbind(ECX_dat, HIP_dat)
  
  p <- plot_DMP(betaMatrix=NULL,phenotypeFile=NULL,position=NULL, interaction=interaction,pathology=pathology,model=model,dat=dat) + facet_grid(~ tissue)
  
  if(!is.null(transcript)){
    p <- p + labs(subtitle = position) + mytheme
    gene_track <- plot_gene_track(ECXbetaMatrix, ECXphenotypeFile, position, model, gene, transcript)
    output <- plot_grid(gene_track,p,nrow=2, rel_heights = c(0.3,0.7))
  }else{
    output <-  p + labs(subtitle = bquote(italic(.(gene)) ~ "(" * .(position) * ")")) + mytheme
  }
  
  return(output)
}


## ------ manhattan plots -----

plot_manhattan <- function(rrbs, array, mode){
  
  if(mode == "Genotype"){
    manhattanPlot <- rrbs %>% dplyr::select(Position, p.val.Genotype) 
    manhattanPlotArray <- array %>% dplyr::select(X,PrZ.GenotypeTG) %>% 
      merge(., mm10_Manifest, by.x = "X", by.y = 0)
  }else{
    manhattanPlot <- rrbs %>% dplyr::select(Position, p.val.Pathology) 
    manhattanPlotArray <- array %>% dplyr::select(X,PrZ.Pathology) %>% 
      merge(., mm10_Manifest, by.x = "X", by.y = 0)
  }
  manhattanPlot$CHR <- stringr::word(manhattanPlot$Position,c(1),sep=stringr::fixed(":"))
  manhattanPlot$BP <- stringr::word(manhattanPlot$Position,c(2),sep=stringr::fixed(":"))
  manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrX")]<-"chr23"
  manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrY")]<-"chr24"
  manhattanPlot$CHR <- as.numeric(str_remove(manhattanPlot$CHR,"chr"))
  manhattanPlot$BP <- as.numeric(manhattanPlot$BP)
  manhattanPlot <- manhattanPlot %>% mutate(Platform = "RRBS")
  #manhattanPlot <- merge(manhattanPlot, rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position") 
  
  manhattanPlotArray$CHR <- stringr::word(manhattanPlotArray$position,c(1),sep=stringr::fixed(":"))
  manhattanPlotArray$BP <- stringr::word(manhattanPlotArray$position,c(2),sep=stringr::fixed(":"))
  manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrX")]<-"chr23"
  manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrY")]<-"chr24"
  manhattanPlotArray$CHR <- as.numeric(str_remove(manhattanPlotArray$CHR,"chr"))
  manhattanPlotArray$BP <- as.numeric(manhattanPlotArray$BP)
  manhattanPlotArray <- manhattanPlotArray  %>% mutate(Platform = "Array")
  
  if(mode == "Genotype"){
    manhattanPlotArray <- manhattanPlotArray %>% dplyr::select(position, PrZ.GenotypeTG, CHR, BP, Platform)
  }else{
    manhattanPlotArray <- manhattanPlotArray %>% dplyr::select(position, PrZ.Pathology, CHR, BP, Platform)
  }
  colnames(manhattanPlotArray) <- colnames(manhattanPlot)
  
  mergedManhattanPlot <- rbind(manhattanPlotArray,manhattanPlot)
  
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
  axisdf <- axisdf %>% filter(CHR != "NA")
  
  donRRBS <- merge(don[don$Platform == "RRBS",], 
                   rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  donArray <- merge(don[don$Platform == "Array",], 
                    rTg4510_array_anno$ECX$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  don <- rbind(donRRBS, donArray)
  if(mode == "Genotype"){
    don <- don %>% arrange(p.val.Genotype)
  }else{
    don <- don %>% arrange(p.val.Pathology)
  }
  top <- unique(don$ChIPseeker_GeneSymbol)[1:200]
  don <- don %>% mutate(label = ifelse(ChIPseeker_GeneSymbol %in% top, ChIPseeker_GeneSymbol, NA))
  don <- don %>% mutate(label_unique = if_else(duplicated(label), NA_character_, label))
  
  don <- don %>% filter(!is.na(CHR))
  don$Platform <- as.factor(don$Platform)
  if(mode == "Genotype"){
    p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Genotype))) 
  }else{
    p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Pathology))) 
  }
  p <- p +
    
    # Show all points
    geom_point(aes(color = Platform, fill = as.factor(Platform)), 
               size = 1.3, shape = 21) +  # Adjust alpha to 0.6 for better visibility
    scale_color_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
    scale_fill_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    #geom_text_repel(data = don[!is.na(don$label_unique), ], 
    #                aes(label = label_unique), 
    #                size = 3, 
    #                box.padding = 0.3, 
    #                point.padding = 0.5) +
    
    # Custom the theme:
    mytheme +
    theme(legend.position="none") +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))))
  
  return(p)
  
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

## ------ Effect size comparisons -----

# all sites
effectSizeComparisons <- function(beta_1, beta_2, platform, model, tissue, animal=NULL){
  
  if(platform == "Array"){
    if(length(colnames(beta_1)["Position" == colnames(beta_1)]) == 0){
      beta_1 <- beta_1 %>% mutate(Position = X)
    }
    if(length(colnames(beta_2)["Position" == colnames(beta_2)]) == 0){
      beta_2 <- beta_2 %>% mutate(Position = X)
    }
  }
  
  if(platform == "Array" & model == "Genotype"){
    cols = c("Position","Betas.GenotypeTG")
  } else if(platform == "RRBS" & model == "Genotype"){
    cols = c("Position","estimate.Genotype")
  } else if(platform == "Array" & model == "Pathology"){
    cols = c("Position","Betas.Pathology")
  } else {
    cols = c("Position","estimate.Pathology")
  }
  
  
  dat <- merge(beta_1[,cols], beta_2[,cols], by = cols[1])
  print(dat)
  
  if(nrow(dat) > 5){
    
    test<- cor.test(dat[[2]], dat[[3]])
    print(test)
    r <- round(test$estimate, 2)
    p <- signif(test$p.value, 3)
    label <- paste0("r = ", r, ", p = ", p)
    
    if(tissue == "HIP"){
      
      if(is.null(animal)){
        print("need to specify whether rTg4510 or J20 for labels")
      }
      
      colnames(dat) <- c("X_Position","ECX", "HIP")
      p <- ggplot(dat, aes(x = ECX, y = HIP)) + geom_point() + 
        theme_classic() +
        labs(x = paste0(animal," ECX"), y = paste(animal, " HIP"), subtitle =  paste0(model,"-associated effect size using ", platform))
      
    }else{
      
      colnames(dat) <- c("X_Position","rTg4510", "J20")
      p <- ggplot(dat, aes(x = rTg4510, y = J20)) + geom_point() + 
        theme_classic() +
        labs(x = "rTg4510", y = "J20", subtitle =  paste0(model,"-associated effect size in ECX using ", platform))
      
    }
    
    p <- p + annotate("text", x = min(dat[[2]]), y =  max(dat[[3]]), label = label, hjust = 0, vjust = 0) 
    return(p)
  }else{
    p <- NULL
    print("Not sufficient common observations")
  }
  
}

