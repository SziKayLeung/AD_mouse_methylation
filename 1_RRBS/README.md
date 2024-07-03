## 1. Bismark  
Bismark alignment using Bowtie2 and methylation extraction. Reference paths etc in ```preprocessing.config```
> [!NOTE]
> Methylation from bismark was further merged for adjacent probes: ```coverage2cytosine --merge_CpG --discordance 100```

## 2. BiSeq smoothing and predicting methylation 
> [!TIP]
> to run using r script. ```predictMeth``` takes a lot of memory and time to run, therefore split into chunks of 100 and merged downstream. 

**Output**  
- raw RRBS: ```0_ZenOutput/1_processed/rrbs/biseq_raw_<model>.RData```  
- smoothed RRBS: ```0_ZenOutput/1_processed/rrbs/<model>_RRBS_SmoothBetas.RData```

## 3. Differential methylation position analysis using BiSeq regression 
> [!TIP]
> to run by submitting bash script that runs r script (genotype, pathology)

**Output**  
- DMP results for
  + Genotype and Interaction: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_betaResultsDMP_GenotypeInteraction.RData```
  + Pathology: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_betaResultsDMP_Pathology.RData```
  + Merged Genotype, Interaction, Pathology as list: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_betaResultsDMPs.RData```

- significant DMP results (FDR < 0.05) for
  + Genotype and Interaction: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_sigResultsDMP_GenotypeInteraction.RData```  
  + Pathology: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_sigResultsDMP_Pathology.RData```
  + Merged Genotype, Interaction, Pathology as list: ```0_ZenOutput/2_differentialAnalysis/rrbs/<model>_sigResultsDMPs.RData```

 ### 4. Annotation using ChIPseeker
 **Output**  
 - annotated significant DMP results (FDR < 0.05) saved as list for Genotype, Interaction, Pathology:
   ```0_ZenOutput/3_annotated/rrbs/<model>_rrbs_annoSigResultsDMPs.RData```
   
Output merges the DMP results with the chipseeker output (also used reference gtf as separate validation of annotation; output have "Ref" prefix)
Column names of annotated genotype dataframe:
```
 [1] "position"                     "Position"                     "FDR_adj_genotype"             "p.val.Genotype"               "meth.group1.WT"              
 [6] "meth.group2.TG"               "meth.diff.Genotype"           "estimate.Genotype"            "std.error.Genotype"           "pseudo.R.sqrt"               
[11] "FDR_adj_age"                  "p.val.Age"                    "estimate.Age"                 "std.error.Age"                "FDR_adj_interaction"         
[16] "p.val.Interaction"            "estimate.Interaction"         "std.error.Interaction"        "p.val.modelLRT"               "ChIPseeker_Annotation"       
[21] "RefGtf_GeneEnsembl"           "RefGtf_GeneSymbol"            "RefGtf_TransEnsembl"          "ChIPseeker_GeneEnsembl"       "ChIPseeker_GeneSymbol"       
[26] "ChIPseeker_GeneName"          "ChIPseeker_TransEnsembl"      "distanceToTSS"                "genic"                        "Intergenic"                  
[31] "Promoter"                     "fiveUTR"                      "threeUTR"                     "Exon"                         "Intron"                      
[36] "downstream"                   "distal_intergenic"            "InGeneBodyOr1500bpTSS"        "InGeneBodyOr1500bpTSS_SYMBOL" "Consensus_GeneSymbol"
```
