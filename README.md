# Methylation profiling of rTg4510 & J20 AD mouse model

This directory is a repository of scripts pertaining to methylation profiling (RRBS, Vertebrate array) of rTg4510 and J20 mice, and form the source of the paper: **Methylomic signatures of tau and beta-amyloid in transgenic mouse models of Alzheimer’s disease neuropathology** by I.Castanho, SK.Leung, ..E.Hannon,J.Mill.

## **Summary**
**Tissue**: 
- rTg4510 Entorhinal Cortex (ECX) and Hippocampus (HIP): n = 31 WT, n = 30 TG, across ages 2, 4, 6, and 8 months
- J20 Entorhinal Cortex (ECX) and Hippocampus (HIP): n = 31 WT, n = 32 TG, ages 6, 8, 10, and 12 months

**Methylation Profiling**:   
- DNA methylation in rTg4510 ECX and J20 ECX was profilied using   
  (i) reduced representation bisulfite sequencing (RRBS), and   
  (ii) DNA methylation microarray
- DNA methylaton in rTg4510 and J20 HIP profiled using only DNA methylation microarray.

**Validation**:
- Pyrosequencing of rTg4510 Ank1 and PrnP differentially methylated sites

## **Bioinformatics pipeline**
> [!NOTE]
> ```import.config``` contains directory path input/output and metadata

1. [RRBS](https://github.com/SziKayLeung/AD_mouse_methylation/tree/dev/1_RRBS) (rTg4510 ECX, J20 ECX)
    + Bismark coverage (note merging of adjacent probes)
    + BiSeq smoothing and predicting methylation across all samples
    + BiSeq differential methylation analysis  
         Genotype: methylation ~ Genotype + Age + Genotype*Age  
         Pathology: methylation ~ Pathology  
    + Annotate DMPs using ChIPseeker  

2. [Array](https://github.com/SziKayLeung/AD_mouse_methylation/tree/dev/2_Array) (rTg4510 ECX, J20 ECX, rTg4510 HIP, J20 HIP)  
   + preprocessing: QC, normlisation  
   + mixed-effects beta regression for differential methylation analysis  
       Genotype: methylation ~ Genotype + Age_months + Genotype*Age_months + Chip_ID  
       Pathology: methylation ~ Pathology + Chip_ID  
   + Annotate DMPs using ChIPseeker
     
3. [Merge](https://github.com/SziKayLeung/AD_mouse_methylation/tree/dev/3_ArrayRRBSComparison) significant DMP (FDR < 0.05) from array and RRBS for rTg4510 ECX and J20 ECX analysis  
