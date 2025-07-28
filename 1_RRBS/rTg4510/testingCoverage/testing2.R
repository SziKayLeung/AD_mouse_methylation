#for i in *bismark.cov*;do 
#  sample=$(basename "$i" | cut -d "_" -f 2 )  
#  echo $sample
#  grep -w 23023240 $i > ${sample}_23023240.cov
#  grep -w 23023241 $i > ${sample}_23023241.cov
#done

library("dplyr")
library("stringr")
library("ggplot2")
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/rTg4510/import.config")
message("Import bismark coverage files after merging CpG")
importSiteCoverage <- function(position){
  files <- list.files(dirnames$bismark, pattern = paste0(position,".cov"), all.files=T, full.names = T)

  inputFiles <- lapply(files, function(x) 
    if(file.size(x) != 0){
      read.table(x, col.names = c("chr","start","end","methylated_perc","methlatedReads","unmethylatedReads"))
    }
    )
  
  names(inputFiles) <-  lapply(files, function(x) if(file.size(x) != 0){word(basename(x),c(1),sep=fixed("_"))})
  merged <- bind_rows(inputFiles, .id = "sample")
  return(merged)
}

Cov23023240 <- importSiteCoverage("23023240")
Cov23023241 <- importSiteCoverage("23023241")

merge(phenotype$rTg4510, rbind(Cov23023240,Cov23023241), by.x = 0, by.y = "sample") %>% 
  ggplot(., aes(x = Genotype, y = methylated_perc)) + facet_grid(~start) +
  geom_boxplot()

manualMerge <- rbind(Cov23023240,Cov23023241) %>% group_by(sample) %>% summarise(methylatedReads = sum(methlatedReads), 
                                                                  unmethylatedReads = sum(unmethylatedReads)) %>% 
  mutate(methylated_perc = methylatedReads/(unmethylatedReads + methylatedReads) * 100) %>%
  merge(., phenotype$rTg4510[,c("Genotype","Age_months")], by.y = 0, by.x = "sample")

SKLMerged <- RRBS_rawbetas_SKL[c("chr8:23023240"),] %>% reshape2::melt(variable.name = "sample", value.name = "methylationSKL")
final <- merge(manualMerge, SKLMerged, by = "sample")


ggplot(manualMerge, aes(x = Genotype, y = methylated_perc)) + 
  geom_boxplot()

merge(phenotype$rTg4510, rbind(Cov23023240,Cov23023241), by.x = 0, by.y = "sample") %>% 
  ggplot(., aes(x = Genotype, y = methylated_perc)) + geom_boxplot()


dat <- merge(tidyr::spread(rbind(Cov23023240) %>% select(sample, methylated_perc),key = sample, value = methylated_perc) %>% 
  reshape2::melt(variable.name = "sample", value.name = "methy23023240"),
  tidyr::spread(rbind(Cov23023241) %>% select(sample, methylated_perc),key = sample, value = methylated_perc) %>% 
    reshape2::melt(variable.name = "sample", value.name = "methy23023241"), all = T) %>%
  mutate(discordant = abs(methy23023240 - methy23023241)) 

res <- cor.test(dat$methy23023240,dat$methy23023241)
ggplot(dat, aes(x = methy23023240, y = methy23023241)) + geom_point() + xlim(0,100) + ylim(0,100) +
  geom_abline() +
  annotate("text", x=15, y=100, label= paste0("r = ", round(res$estimate,4))) + 
  labs(x = "Methylation; chr8: 23023240", y = "Methylation; chr8: 23023241")

res$estimate

final <- merge(final, dat[,c("sample","discordant")], by = "sample")  %>% arrange(discordant)                                                              
