#!/bin/sh
#PBS -V #export all enviroment variables to the batch job
#PBS -q pq #submit to the serial queue
#PBS -l walltime=3:00:00 # maximum wall time for the job
#PBS -A Research_Project-MRC190311 #research project to submit under
#PBS -l nodes=1:ppn=16 #Number of processors 
#PBS -m e -M a.n.dahir@exeter.ac.uk # email me at job completion 

module load R/3.6.0-foss-2019a

cd /gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/Array/ # Setting Working directory to where script is located 

Rscript "/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/Array/MixedModelsEffectBetaReg_rTg4510.R"  ## execute R script

Rscript "/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/Array/MixedModelsEffectBetaReg_J20.R"

