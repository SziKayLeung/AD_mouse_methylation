## SupplementaryFigures
## correlations between rTg4510 vs J20 all sites and sig plots
## correlations between HIP and ECX
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config.R"))

# rTg4510 vs J20 
message("Number of common significant sites between rTg4510 and J20: ", 
        length(intersect(sigRes$rTg4510$Genotype[sigRes$rTg4510$Genotype$Platform == "Array",c("Position")],
          sigRes$J20$Genotype[sigRes$J20$Genotype$Platform == "Array",c("Position")])))

message("Number of common significant sites between rTg4510 and J20: ", 
        length(intersect(sigRes$rTg4510$Pathology[sigRes$rTg4510$Pathology$Platform == "Array",c("Position")],
                         sigRes$J20$Pathology[sigRes$J20$Pathology$Platform == "Array",c("Position")])))

merge(sigRes$rTg4510$Genotype,
      sigRes$J20$Genotype, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point()

merge(sigRes$rTg4510$Pathology,
      sigRes$J20$Pathology, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point()


# hippocampus significant results
p1 <- merge(sigResArrayECX$rTg4510$Genotype,
      sigResArrayHIP$rTg4510$Genotype, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point() +
  labs(x = "ECX rTg4510 effect size  - genotype", y = "HIP rTg4510 effect size - genotype")
  
p2 <- merge(sigResArrayECX$J20$Genotype,
            sigResArrayHIP$J20$Genotype, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point() +
  labs(x = "ECX J20 effect size - genotype", y = "HIP J20 effect size - genotype")

p3 <- merge(sigResArrayECX$rTg4510$Pathology,
            sigResArrayHIP$rTg4510$Pathology, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point() +
  labs(x = "ECX rTg4510 effect size - pathology", y = "HIP rTg4510 effect size - pathology")

p4 <- merge(sigResArrayECX$J20$Pathology,
            sigResArrayHIP$J20$Pathology, by = "Position") %>% 
  ggplot(., aes(x = delta.x, y = delta.y)) + geom_point() +
  labs(x = "ECX J20 effect size - pathology", y = "HIP J20 effect size - pathology")


## ------------------------- all sites

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

# rTg4510 vs J20
p1 <- effectSizeComparisons(rTg4510_array_results$Genotype, J20_array_results$Genotype, "Array", "Genotype", "ECX")
p2 <- effectSizeComparisons(rTg4510_array_results$Pathology, J20_array_results$Pathology, "Array", "Pathology", "ECX")
p3 <- effectSizeComparisons(rTg4510_rrbs_results$Genotype, J20_rrbs_results$Genotype, "RRBS", "Genotype","ECX")
p4 <- effectSizeComparisons(rTg4510_rrbs_results$Pathology, J20_rrbs_results$Pathology, "RRBS", "Pathology","ECX")

pdf(paste0(output,"Figures/rTg4510vsJ20.pdf"), width = 10, height = 15)
plot_grid(p1,p2,p3,p4,labels = c("A","B","C","D"))
dev.off()

# significant 
effectSizeComparisons(rTg4510_array_sig$ECX$Genotype, J20_array_sig$ECX$Genotype, "Array", "Genotype", "ECX")
effectSizeComparisons(rTg4510_array_sig$ECX$Pathology, J20_array_sig$ECX$Pathology, "Array", "Pathology", "ECX")

# HIP vs ECX
p5 <- effectSizeComparisons(rTg4510_array_results$Genotype, rTg4510_HIP_array_results$Genotype, "Array", "Genotype", "HIP", animal="rTg4510")
p6 <- effectSizeComparisons(rTg4510_array_results$Pathology, rTg4510_HIP_array_results$Pathology, "Array", "Pathology", "HIP", animal="rTg4510")
p7 <- effectSizeComparisons(rTg4510_array_sig$ECX$Genotype, rTg4510_array_sig$HIP$Genotype, "Array", "Genotype", "HIP", animal = "rTg4510")
p8 <- effectSizeComparisons(rTg4510_array_sig$ECX$Pathology, rTg4510_array_sig$HIP$Pathology, "Array", "Pathology", "HIP", animal = "rTg4510")
pdf(paste0(output,"Figures/rTg4510ECXvsHIP.pdf"),  width = 10, height = 15)
plot_grid(p5,p6,p7,p8,labels = c("A","B","C","D"))
dev.off()


p9 <- effectSizeComparisons(J20_array_results$Genotype, J20_HIP_array_results$Genotype, "Array", "Genotype", "HIP", animal="J20")
p10 <- effectSizeComparisons(J20_array_results$Pathology, J20_HIP_array_results$Pathology, "Array", "Pathology", "HIP", animal="J20")
p11 <- effectSizeComparisons(J20_array_sig$ECX$Genotype, J20_array_sig$HIP$Genotype, "Array", "Genotype", "HIP", animal = "J20")
p12 <- effectSizeComparisons(J20_array_sig$ECX$Pathology, J20_array_sig$HIP$Pathology, "Array", "Pathology", "HIP", animal = "J20")
pdf(paste0(output,"Figures/J20ECXvsHIP.pdf"),  width = 10, height = 15)
plot_grid(p9,p10,p11,p12,labels = c("A","B","C","D"))
dev.off()

dat <- merge(rTg4510_array_results$Genotype[, c("X", "Betas.GenotypeTG")], 
             rTg4510_HIP_array_results$Genotype[, c("X", "Betas.GenotypeTG")], by = "X") %>% 
  mutate(sig_ECX = ifelse(X %in% rTg4510_array_sig$ECX$Genotype$cpg, TRUE, FALSE)) %>% 
  mutate(sig_HIP = ifelse(X %in% rTg4510_array_sig$HIP$Genotype$cpg, TRUE, FALSE)) %>%
  mutate(sig = case_when(
    sig_ECX & sig_HIP ~ "both",
    sig_ECX ~ "ECX",
    sig_HIP ~ "HIP",
    TRUE ~ "none"
  ))

ggplot(dat, aes(x = Betas.GenotypeTG.x, y = Betas.GenotypeTG.y, colour = sig)) + geom_point() +
  labs(x = "ECX", y = "HIP")
