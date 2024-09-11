# rTg4510
total = 20023
up = 15558
down = 4465

# J20
total = 16473
up = 11672
down = 4801










# plot signif genes from Tg4510 genotype effect againts same genes in J20
threshold <- 0.05

#ntotal <- sum(Tg4510_genotype$padj < threshold, na.rm=TRUE) # number of sig genes I previously identified? should be equal or less...
ntotal <- total # number of sig genes I previously identified? should be equal or less...
#nconsistent <- sum(sign(Tg4510_genotype$log2FoldChange[which(Tg4510_genotype$padj < threshold)]) == sign(J20_genotype$log2FoldChange[which(Tg4510_genotype$padj < threshold)]), na.rm=TRUE) # list of T and F - sum will give me the number of times T
nconsistent <- up # list of T and F - sum will give me the number of times T

test <- binom.test(nconsistent, ntotal, 0.5)

pvalue <- test[3]
pvalue <- signif(as.numeric(pvalue), 3)
pvalue # p-value = 0.468
probability <- test[5]
probability <- signif(as.numeric(probability), 3)
probability # p = 0.468