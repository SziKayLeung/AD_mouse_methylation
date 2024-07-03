## 1. QC and normalisation 
*Processed by A.Dahir*
- checked no batch effects
- removed low intensity samples
- removed samples with low bisulphite conversion
- normalisation using ```HorvathMammal40.CanonicalManifest.3.2019.sesame.csv```

**Output**  
Normalised rTg4510 and J20 ECX and HIP: ```0_ZenOutput/1_processed/array/Normalised_Data_Sesame.rdat```  
Normalised rTg4510 ECX (Subsetted): ```0_ZenOutput/1_processed/array/Normalised_rTg4510_array_ECX.RData```  
Normalised J20 ECX (Subsetted): ```0_ZenOutput/1_processed/array/Normalised_J20_array_ECX.RData```  
  
## 2. Differential methylation position analysis using mixed-model beta regression 
*Processed by E.Walker*
> [!Note]
> Analysis were performed by submitting bash script that runs r script  
> ECX, HIP: genotype, interaction: ```methylation ~ Genotype + Age_months + Genotype*Age_months + Chip_ID```  
> ECX, HIP: pathology: ```methylation ~ Pathology + Chip_ID```    
> Tissue specific: ```methylation ~ Tissue + Tissue*Genotype + Genotype + Age_months + (1|MouseID)```  
> For Zenodo upload, output were moved and saved under different name to avoid confusion (check ```resaving_files.R```)

**Output**  
- DMP results for
  + ECX (List of genotype, interaction, pathology results): ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_ECX_allResultsDMPs.RData```
  + HIP (List of genotype, interaction, pathology results): ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_HIP_allResultsDMPs.RData```
  + Tissue difference: ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_tissue_allResultsDMPs.RData```

- significant DMP results (FDR < 0.05) for
  + ECX (List of genotype, interaction, pathology results): ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_ECX_sigResultsDMPs.RData```
  + HIP (List of genotype, interaction, pathology results): ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_HIP_sigResultsDMPs.RData```
  + Tissue difference: ```0_ZenOutput/2_differentialAnalysis/array/<model>_array_tissue_sigResultsDMPs.RData```
    
DMP results were merged with CpG from manifest, created co-ordinate column, FDR adjusted using script ```functions/finalise_betaRegressionArray.R```.   
Output therefore complementary to RRBS output from DMP analysis

 ## 3. Annotation using ChIPseeker
 **Output**  
 - annotated significant DMP results (FDR < 0.05) saved as list for Genotype, Interaction, Pathology:
   ```0_ZenOutput/3_annotated/array/<model>_array_annoSigResultsDMPs.RData"```
   
Output merges the DMP results with the chipseeker output
Column names of annotated genotype dataframe:
```
 [1] "position"                     "cpg"                          "Bp"                           "Gene_Symbol"                  "Betas.Baseline"              
 [6] "SE.Baseline"                  "Z.Baseline"                   "PrZ.Baseline"                 "Betas.GenotypeTG"             "SE.GenotypeTG"               
[11] "Z.GenotypeTG"                 "PrZ.GenotypeTG"               "Betas.Age_months"             "SE.Age_months"                "Z.Age_months"                
[16] "PrZ.Age_months"               "Betas.GenotypeTG.Age_months"  "SE.GenotypeTG.Age_months"     "Z.GenotypeTG.Age_months"      "PrZ.GenotypeTG.Age_months"   
[21] "Betas.Chip_ID204027420024"    "SE.Chip_ID204027420024"       "Z.Chip_ID204027420024"        "PrZ.Chip_ID204027420024"      "Betas.Chip_ID204027420027"   
[26] "SE.Chip_ID204027420027"       "Z.Chip_ID204027420027"        "PrZ.Chip_ID204027420027"      "Betas.Chip_ID204027420037"    "SE.Chip_ID204027420037"      
[31] "Z.Chip_ID204027420037"        "PrZ.Chip_ID204027420037"      "Position"                     "FDR_adj_GenotypeAge"          "FDR_adj_Genotype"            
[36] "seqnames"                     "start"                        "end"                          "width"                        "strand"                      
[41] "ChIPseeker_Annotation"        "geneChr"                      "geneStart"                    "geneEnd"                      "geneLength"                  
[46] "geneStrand"                   "geneId"                       "ChIPseeker_TransEnsembl"      "distanceToTSS"                "ChIPseeker_GeneEnsembl"      
[51] "ChIPseeker_GeneSymbol"        "ChIPseeker_GeneName"          "genic"                        "Intergenic"                   "Promoter"                    
[56] "fiveUTR"                      "threeUTR"                     "Exon"                         "Intron"                       "downstream"                  
[61] "distal_intergenic"            "InGeneBodyOr1500bpTSS"        "InGeneBodyOr1500bpTSS_SYMBOL"
```
