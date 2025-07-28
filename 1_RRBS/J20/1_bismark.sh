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
#SBATCH --array=0-62%10 # 63 samples
#SBATCH --output=1_bismark-%A_%a.o.o
#SBATCH --error=1_bismark-%A_%a.o.e

# 18.06.2024: J20 samples (n=63)

##-------------------------------------------------------------------------

# load 
module load Miniconda2
source activate nanopore
module load Bowtie2
source ./preprocessing.config 

# single sequencing reads
f1=($(ls ${TRIMDIR}/${sampleName}*_trimmed*.fq.gz*))
sample=${f1[${SLURM_ARRAY_TASK_ID}]} 
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
gunzip *bismark.cov.gz*

# alignment, already aligned from Thea (backed up on hard-drive)
echo 'Merge CpG all independent of discordance rate'
coverage2cytosine ${sampleName}.bismark.cov --merge_CpG --discordance 100 --output ${sampleName} --genome_folder ${REFPATH} &>> ${sampleName}_coverage2cytosine.log 