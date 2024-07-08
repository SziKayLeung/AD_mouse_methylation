#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --array=0-6 # 7 samples
#SBATCH --output=1b_remaining_bismark-%A_%a.o.o
#SBATCH --error=1b_remaining_bismark-%A_%a.o.e

# 13.05.2024: run bismark_genome_preparation 
# 16.05.2024: run bismark on one sample (1.5hr)
# https://github.com/FelixKrueger/Bismark/issues/669 
# 20.05.2024: rTg4510: 54 samples
# 04.06.2024: running with 100% discordance at coverage2cytosine i.e. merge everything
# 10.06.2024: run coverage2cytosine for missing 6 samples from trimmed directory (using coverage files backed up from drive)

##-------------------------------------------------------------------------

# load 
module load Miniconda2
source activate nanopore
module load Bowtie2
source ./preprocessing.config 

# single sequencing reads
f1=($(ls ${ALIGNDIR}/*cov.gz*))
sample=${f1[${SLURM_ARRAY_TASK_ID}]} 
sampleName=$(basename "$sample" | cut -d "." -f 1 )

##-------------------------------------------------------------------------

cd ${ALIGNDIR}

# bismark, alignment with bowtie2
echo 'Merge CpG all independent of discordance rate'
coverage2cytosine ${sampleName}.bismark.cov.gz --merge_CpG --discordance 100 --output ${sampleName} --genome_folder ${REFPATH} &>> ${sampleName}_coverage2cytosine.log 