# Methylation profiling of rTg4510 & J20 AD mouse model

This directory is a repository of scripts pertaining to methylation profiling (RRBS, Vertebrate array) of rTg4510 and J20 mice, and form the source of the paper: **Methylomic signatures of tau and amyloid-beta in transgenic mouse models of Alzheimer’s disease neuropathology** by SK.Leung,E.M.Walker...I.Castanho,J.Mill.

## **Summary**
**Tissue**: 
- rTg4510 Entorhinal Cortex (ECX) and Hippocampus (HIP): n = 31 WT, n = 30 TG, across ages 2, 4, 6, and 8 months
- J20 Entorhinal Cortex (ECX) and Hippocampus (HIP): n = 32 WT, n = 31 TG, ages 6, 8, 10, and 12 months

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

## **Manuscript**

To recreate main figures, supplementary figures and tables, the processed data from `0_PaperFigs/paper_import.config.R` needs to be imported.  
Note: the paths to this processed data is currently directed to S.Leung's local drive, therefore the following paths need to be updated:
```
# this github repo
scriptDir = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/AD_mouse_methylation/"
# processed data 
output = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/PaperOutput/"
# github repo: https://github.com/SziKayLeung/LOGen
LOGEN_ROOT = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/LOGen/"
```
The `scriptDir` variable in the MainFigures.R, SupplementaryFigures.R, Stats.R, Tables.R also needs to be updated to the location of this github repo in order to source `paper_import.config.R`.

The data from annotation onwards can be downloaded from Zenodo (https://zenodo.org/records/15741354). If processed data is required, please email S.K.Leung@exeter.ac.uk for access. 
All data is stored on ISCA: `/lustre/projects/Research_Project-MRC148213/lsl693/rrbs_ad_mice/0_ZenOutput`
