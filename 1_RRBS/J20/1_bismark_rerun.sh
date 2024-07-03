#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=30:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --array=3-16%5 # all 17 samples
#SBATCH --output=1_bismark_rerun-%A_%a.o.o
#SBATCH --error=1_bismark_rerun-%A_%a.o.e

# 18.06.2024: J20 samples (n=63); error with one of the sample therefore re-run 2605_A17_r1

##-------------------------------------------------------------------------

# load 
module load Miniconda2
source activate nanopore
module load Bowtie2
source ./preprocessing.config 

# single sequencing reads
reRunSamples=(A18 A23 F19 G24 H18 H22 H23 H24 I20 J17 J19 J20 J21 J22 J23 J24 A17)
reRunSample=${reRunSamples[${SLURM_ARRAY_TASK_ID}]}
sample=${DATADIR}/2_trimmed/2605_${reRunSample}_r1_trimmed.fq.gz
sampleName=$(basename "$sample" | cut -d "." -f 1 )


##-------------------------------------------------------------------------

# bismark, alignment with bowtie2
# already aligned from processing with rTg4510
if [ -d "${REFPATH}/Bisulfite_Genome" ]; then
echo "Bismark genome preparation performed"
else
  echo "Running bismark genome preparation"
bismark_genome_preparation --path_to_aligner ${BOWTIE2} ${REFPATH} &>> ${ALIGNDIR}/bismark_genome_prep.log
fi

# alignment
cd ${ALIGNDIR}
echo 'run Bismark on' ${sampleName}
bismark --genome ${REFPATH} ${sample} --basename ${sampleName} -o ${ALIGNDIR} &>> ${ALIGNDIR}/${sampleName}_bismark.log

echo 'Extracting methylation'
# single end
bismark_methylation_extractor --single-end ${sampleName}.bam --bedGraph -o ${ALIGNDIR} --merge_non_CpG --comprehensive &>> ${ALIGNDIR}/${sampleName}_bismark_extractor.log 

# alignment, already aligned from Thea (backed up on hard-drive)
echo 'Merge CpG all independent of discordance rate'
coverage2cytosine ${sampleName}.bismark.cov.gz --merge_CpG --discordance 100 --output ${sampleName} --genome_folder ${REFPATH} &>> ${sampleName}_coverage2cytosine.log 
