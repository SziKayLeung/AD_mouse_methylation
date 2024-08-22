#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../output/5_calculate_vario.o
#SBATCH --error=../output/5_calculate_vario.e

module load R/4.2.2-foss-2022b
#Rscript /lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/1_RRBS/rTg4510/5_identify_DMR.R
Rscript /lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/1_RRBS/rTg4510/5_calculate_vario.R