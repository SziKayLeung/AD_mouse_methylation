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
#SBATCH --array=0-7 # 8 rates
#SBATCH --output=1a_bismark-%A_%a.o.o
#SBATCH --error=1a_bismark-%A_%a.o.e

# 30/05/2024: testing different discordant rates on 4 samples

##-------------------------------------------------------------------------

# load 
module load Miniconda2
source activate nanopore
module load Bowtie2
source ./preprocessing.config 

# single sequencing reads
f1=($(ls ${TRIMDIR}/${sampleName}*_trimmed*.f*))
discordanceRates=(30 40 50 60 70 80 90 100)
discordanceRate=${discordanceRates[SLURM_ARRAY_TASK_ID]}

echo "Testing out: $discordanceRate%"

# trialling 4 samples
cd ${ALIGNDIR}/sensitivity
for i in {1..4}; do 

  echo $i
  sample=${f1[$i]} 
  sampleName=$(basename "$sample" | cut -d "." -f 1 )
  echo $sample
  coverage2cytosine ${ALIGNDIR}/${sampleName}.bismark.cov --merge_CpG --discordance $discordanceRate --output ${sampleName}_${discordanceRate} --genome_folder ${REFPATH} &>> ${sampleName}_coverage2cytosine_${discordanceRate}.log 

done

