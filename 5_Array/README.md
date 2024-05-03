Array

QC : 1_QC_Mouse_Array.Rmd 
	 - Runs quality control on raw data and produces normalised rdat. Removed failed samples.

Array Summary: 2_Array_summary.Rmd
	 - Summaries the reader on what about the custom array, types of probes (unique vs ambiguous), other species it maps into, overlap with EPIC/450K array.
	 - Also does the gene annotation on the array and compares number of site that map onto the same gene

Data Exploration: 3_Explore_Array_Data.Rmd
	 - Explore the data e.g PCA plots

Mixed Models: 
4_MEM_XX_Tissue_sbatch.sh; 4a_MixedModelsEffectBetaReg_XX_Tissue_EW.R
   - Run separately on rTg4510 and J20 
   - model: ~ Tissue + Tissue*Genotype + Genotype + Age_months + (1|MouseID)
   - job submission script and Rscript (EW submitted in personal folder)

4b_MixedModelsEffectBetaReg_XX_ECX_EW.R
   - Run separately on rTg4510 and J20 
   - model: ~ Genotype + Age_months + Genotype*Age_months + Chip_ID

Mixed Models Analysis: 5_MEM_analysis; 5a_MEM_analysisFunctions.R
   - Annotate array probes using ChipSeeker
   - Plot output from mixed models regression 

---Archived---
Annotation: mm10_Annotations_chipseeker.R
	 - Script for gene annotation although this is done in the Array_summary.Rmd script
	 - Don't need to run this script.

EWAS: rTg4510_Array_BetaRegression.Rmd + rTg4510_Hip_Array_BetaRegression.Rmd + J20_Array_BetaRegression.Rmd + J20_Hip_Array_BetaRegression.Rmd
	 - Rmd runs through seperate single models of the analysis as well as full model (Genotype + Age + Pathology)

Significant Site Threshold: Array_FDRvsBonf.R
	 - Plots a table of signifcant hits using either P value correction method

Compare Tissues : CompareRawMethylationTissues.R  + Genotype_CompareTissues.R + Interaction_CompareTissues.R + Pathology_CompareTissues.R + InteractionPathology_Comparison.R
	 - Compared the DMPs from single models to the tissues in both the Cortex and Hippocampus

Compare pathology Models: rTg4510_Pathology.R + J20_Pathology.R + PathologyModelsResults.R
	 - These scripts were to check that the pathology term alone differs to the full model 
	

Mixed Models: rTg4510_MixedModelsEffect.Rmd + J20_MixedModelsEffect.Rmd + MixedModelsEffectBetaReg_J20.R + MixedModelsEffectBetaReg_rTg4510.R + MixedModelsEffectPlot_J20.R + MixedModelsEffectPlotsrTg.R
	 - MixedModels are run in the MixedModelsEffectBetaReg* scripts then can be plotted in the MixedModelsEffectPlot*. 
	 - The Rmd summarises the plot R scripts and saves the outputs

Clocks: DNAmAge-Clocks.Rmd + Clock.R + ClockCortex.R + ClockCortexFullModel.R
	 - Rmd summaries/repeats what is done in the individual R scripts. 
	 - Clock.R plots the pdf to visualise age vs predicted age - Also runs hippocampus samples using brain clock and cortex samples using cortex clock
	 - ClockCortex.R runs models using clock for cortex. ClockCortexFullModel.R runs mixed models (combining tissues only)

Functions: BetaRegressionArray.R
	 - Functions to run beta aggression in single and mixed models.