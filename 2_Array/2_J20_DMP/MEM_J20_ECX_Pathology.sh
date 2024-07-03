#!/bin/sh

#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=3:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END
#SBATCH --mail-user=E.M.Walker@exeter.ac.uk # email me at job completion
#SBATCH --job-name=J20_ECX_Path

module load Miniconda2 

source activate EWAS_parallel  # activating environment

cd /gpfs/ts0/projects/Research_Project-191406/EmmaW/Array # Setting Working directory to where script is located 

Rscript "/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/MixedModelsEffectBetaReg_J20_ECX_Pathology_EW.R"  ## execute R script

