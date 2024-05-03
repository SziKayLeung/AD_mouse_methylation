#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/code/RRBS-GitLab/ # set working directory to directory where script is
#SBATCH -p mrcq # submit to the mrc queue
#SBATCH --time=24:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=I.S.Castanho@exeter.ac.uk # email address

## execute R script

module load R/3.6.0-foss-2019a
# module load R/4.0.0-foss-2020a

Rscript 04-J20_RRBS_DMPs-Pathology.R > 04-J20_RRBS_DMPs-Pathology.R-log.r