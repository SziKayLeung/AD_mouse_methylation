
library("tidyr")
library("ggrepel")
library("stringr")


label_colour <- function(var){
  if(var == "Tg4510"){colour = "#00AEC9"}else{
    if(var == "J20"){colour = "#FF5A62"}else{
    }}
  return(colour)
}


##--- Manhattan plots -------

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

merge_manhattan <- function(df1,df2,name1,name2){
  merged <- merge(df1  %>% dplyr::select(Location, mFDR, ChIPseeker_GeneSymbol),
                  df2  %>% dplyr::select(Location, mFDR, SYMBOL),by = "Location",all=T)%>%
    dplyr::select(Location, mFDR.x,mFDR.y,ChIPseeker_GeneSymbol,SYMBOL) %>%
    `colnames<-`(c("Location", name1,name2,"gene1","gene2")) %>%
    mutate(gene1 = as.character(gene1),
           gene2 = as.character(gene2),
           gene = ifelse(!is.na(gene1),gene1,gene2))
  
  merged_prepped <- prepare_manhattan(merged,"none","merge")
  
  return(merged_prepped)
}




plot_manhattan <- function(df, p_col, ylab, title){
  
  # Set p-value trehsholds
  threshold <- 0.05
  sig2 <- 0.01
  
  
  
  p <- manhattan(df,
                 p = p_col,
                 ylim = c(0,6),
                 ylab = ylab,
                 cex = 2, # point size
                 cex.axis = 2.5,
                 cex.lab = 2,
                 col = c("black", color_Tg4510_TG),
                 suggestiveline = -log10(threshold),
                 genomewideline = -log10(sig2),
                 chrlabs = c(1:19, "X"),
                 # annotatePval = 1.947822e-34,
                 las = 2, # rotate x axis labels so they all fit on the plot
                 annotateTop = FALSE,
                 main = title)
  
  return(p)
  
}

plot_manhattan_2 <- function(prep_mhat_df,mouse_model,type,title){
  
  df <- prep_mhat_df %>% dplyr::select(SNP,CHR,BP,contains("mFDR"))
  
  if(type == "single"){
    p <- CMplot(df,plot.type="m",LOG10=TRUE,multracks=FALSE,threshold=c(1e-2),file="jpg",memo="",dpi=300,
                file.output=FALSE,verbose=TRUE,width=14,height=6,
                amplify=TRUE,bin.size=1e6,main=title,
                col=c("black",label_colour(mouse_model)))
    
  }else{
    p <- CMplot(df,plot.type="m",multracks=TRUE,LOG10=TRUE,threshold=c(1e-6,1e-4),file="jpg",memo="",dpi=300,
                file.output=FALSE,verbose=TRUE,width=14,height=6,
                col=c(alpha(label_colour(mouse_model),0.2),label_colour(mouse_model))) +
      print(str(p))
  }
  
  return(p)
}




cumulative <- function(df){
  
  data_cum <- df %>% 
    group_by(CHR) %>% 
    summarise(max_bp = max(BP)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    dplyr::select(CHR, bp_add)
  
  final_cum <- df %>% 
    inner_join(data_cum, by = "CHR") %>% 
    mutate(bp_cum = BP + bp_add)
  
  return(final_cum)
}

plot_manhattan_3 <- function(input_df, tissue, mouse_model,filter){
  # adapted from https://github.com/YinLiLin/CMplot/blob/master/R/CMplot.r
  
  
  if(tissue %in% c("ECX","ECX_upwards")){
    input <- cumulative(gather(input_df, platform, mFDR, mFDR_RRBS:mFDR_Array, factor_key=TRUE)) %>% 
      mutate(anno = ifelse(mFDR < 1e-50, as.character(gene),""))
    input$platform <- str_replace(input$Platform, "mFDR_Array", "Array")
    input$platform <- str_replace(input$Platform, "mFDR_RRBS", "RRBS")
    input$platorm <- as.factor(input$Platform)
    if(filter == "yes"){ input <- input %>% filter(mFDR > 1e-50)}
    
    plot <- ggplot(input, aes(x = bp_cum, y = -log10(mFDR), color = platform,label = anno)) +
      geom_point(aes(shape=platform), alpha = 0.75) +
      scale_shape_manual(values=c(1, 19))
    
  }else{
    input <- input_df %>% mutate(anno = ifelse(mFDR < 1e-50, as.character(gene),""))
    
    if(filter == "yes"){input <- input %>% filter(mFDR > 1e-50)}
    
    plot <- ggplot(input, aes(x = bp_cum, y = -log10(mFDR), color = as.factor(CHR),label = anno)) +
      geom_point(alpha = 0.75)
  }
  
  axis_set <- input %>% 
    group_by(CHR) %>% 
    summarize(center = mean(bp_cum))
  
  axis_edge <- data.frame(input %>% group_by(CHR) %>% summarise(edge = max(bp_cum)))
  
  
  ylim <- abs(floor(log10(min(input$mFDR, na.rm = T)))) + 2
  
  sig <- 5e-2
  gwas_sig <- 1e-3
  
  plot <- plot +  
    geom_text_repel(show.legend = F) + 
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
    geom_hline(yintercept = -log10(gwas_sig), color = "blue", linetype = "dashed") +
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
    ) + guides(shape="none") 
  
  if(tissue == "HIP"){
    if(mouse_model == "Tg4510"){
      altcol = "darkblue"
    }else{
      altcol = "darkred"
    }
    plot <- plot + 
      scale_color_manual(values = rep(c(altcol,label_colour(mouse_model)), 
                                      unique(length(axis_set$CHR)))) + theme(legend.position = "none",
                                                                             axis.text.x=element_blank(),
                                                                             axis.ticks.x=element_blank(),
                                                                             axis.line.x=element_blank())
    
  }
  
  if(tissue == "ECX"){
    plot <- plot + scale_colour_manual(values = c(wes_palette("Zissou1")[3],alpha(label_colour(mouse_model),0.5)), name = "Platform") 
  }
  
  if(tissue %in% c("upwards","ECX_upwards")){
    plot <- plot + scale_color_manual(values = rep(c(label_colour(mouse_model),"black"), 
                                            unique(length(axis_set$CHR))))# + theme(legend.position = "none")
  }
  

  
  return(plot)
  
}
