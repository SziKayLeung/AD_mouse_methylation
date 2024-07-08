#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=80:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../output/3_run_DMP_Pathology.o
#SBATCH --error=../output/3_run_DMP_Pathology.e

module load R/4.2.2-foss-2022b
scriptDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/1_RRBS/J20/
Rscript ${scriptDir}/3b_source_DMP_Pathology.R