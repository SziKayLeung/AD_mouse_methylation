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


