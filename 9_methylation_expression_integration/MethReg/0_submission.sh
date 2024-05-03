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
#SBATCH --array=0-19%5 # 20 chromosomes, 1 - 19, and chr X
#SBATCH --output=MethRegOutput/methRegPromoterTg4510-%A_%a.o
#SBATCH --error=MethRegOutput/methRegPromoterTg4510-%A_%a.e

AllChr=(chr{1..19} chrX)
Chr=${AllChr[${SLURM_ARRAY_TASK_ID}]}

module load R
scriptDir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/9_methylation_expression_integration
#Rscript ${scriptDir}/run_methreg_unsupervised.R -c ${Chr} -t rrbs
Rscript ${scriptDir}/run_methreg_unsupervised.R -c ${Chr} -t array -n arrayUnsupervised