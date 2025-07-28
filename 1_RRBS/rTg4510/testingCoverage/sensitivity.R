
suppressMessages(library("BiSeq"))
suppressMessages(library("stringr"))


sensitivity = list.files(path = "/lustre/projects/Research_Project-MRC148213/lsl693/rrbs_ad_mice/RRBS/rTg4510/3_bismark/sensitivity", 
                         pattern = "CpG_report.merged_CpG_evidence.cov", all.files = T, full.names  = T)


samples <- word(unlist(lapply(sensitivity, function(x) basename(x))),c(2), sep = fixed("_"))
discordance <- word(unlist(lapply(sensitivity, function(x) basename(x))),c(5), sep = fixed("_"))

sensitivity <- lapply(sensitivity, function(x) data.table::fread(x, data.table = F))
names(sensitivity) <- paste0(samples,".",discordance)
sensitivityAdd <- bind_rows(sensitivity, .id = "samples")
sensitivityAdd %>% group_by(samples) %>% tally() %>% 
  mutate(coverage = word(samples,c(2), sep = fixed("."))) %>% 
  mutate(sample = word(samples,c(1), sep = fixed("."))) %>%
  as.data.frame() %>%
  ggplot(., aes(x = as.numeric(coverage), y = n, group = sample)) + geom_point() + geom_line()

dat <- sensitivityAdd %>% group_by(samples) %>% tally() %>% 
  mutate(coverage = word(samples,c(2), sep = fixed("."))) %>% 
  mutate(sample = word(samples,c(1), sep = fixed("."))) %>%
  as.data.frame() 

merge(dat %>% filter(coverage == 100) %>% select(sample,n), dat, by = "sample") %>% 
  mutate(perc_retained = (n.y/n.x) * 100) %>% 
  ggplot(., aes(x = as.numeric(coverage), y = perc_retained, group = sample)) + geom_point() + geom_line() +
  labs(x = "discordance rate", y = "percentage of sites retained")
