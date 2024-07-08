#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../output/3_run_DMP_J20.o
#SBATCH --error=../output/3_run_DMP_J20.e

module load R/4.2.2-foss-2022b
Rscript /lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/1_RRBS/J20/3a_source_DMP_Genotype.R