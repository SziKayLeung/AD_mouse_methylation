#download mouse genome  from EMSEMBLE
cd /mnt/data1/reference_files/mouse_GRCm38.p5
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.7_GRCm38.p5/GCA_000001635.7_GRCm38.p5_genomic.fna.gz 

#latest minor release 07/03/2017 p5

#rename to fa.gz as .fna won't work in bismark


cp GCA_000001635.7_GRCm38.p5_genomic.fna.gz  GCA_000001635.7_GRCm38.p5_genomic.fa.gz 
#genome prepararion with Bismark:

#unzip it 
gzip -d /mnt/data1/reference_files/mouse_GRCm38.p5/GCA_000001635.7_GRCm38.p5_genomic.fa.gz 

#append the spike in sequences to the reference genomes/all/GCA/000/001/635/GCA_000001635

cd /mnt/data1/reference_files/mouse_GRCm38.p5
#append the spike-in sequences to the ref genome fasta file 
 cat Spike_in_Diagenode.fa GCA_000001635.7_GRCm38.p5_genomic.fa > GCA_000001635.7_GRCm38.p5_genomic.spikein.fa
wc GCA_000001635.7_GRCm38.p5_genomic.spikein.fa
wc GCA_000001635.7_GRCm38.p5_genomic.fa
wc Spike_in_Diagenode.fa
#the result should be a sum

#move GCA_000001635.7_GRCm38.p5_genomic.spikein.fa to folder /mnt/data1/reference_files/mouse_GRCm38.p5/appended/


#run bismark on ref genome + spike-in
cd /mnt/data1/Joana/MATRICS/RRBS_test
/mnt/data1/programs/Bismark/bismark_genome_preparation --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1 --verbose /mnt/data1/reference_files/mouse_GRCm38.p5/appended/ 
 


#Use trim_galore and runs fastqc on output - use mobaxterm, fastqc won't open on putty
#remember that the adapters are added before bisulfite conversion, so any C's need to be T's now
#change error rate to 0.2 to allow at least 1 mismatch
#balb13
cd /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A17/
gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A17/raw_illumina_reads/A17_S1_R1_001.fastq.gz

#in this run trim the Diagenode adapter 
trim_galore -q 20  --phred33 --rrbs --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A17/raw_illumina_reads/A17_S1_R1_001.fastq

#in this run trim the Diagenode adapter error rate 20%
trim_galore -a gtagag -q 20 -e 0.2 --phred33 --rrbs --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC2/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore2  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A17/raw_illumina_reads/A17_S1_R1_001.fastq

#balb19

cd /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A24/
gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A24/raw_illumina_reads/A24_S2_R1_001.fastq.gz

#in this run trim the Diagenode adapter
trim_galore -q 20 --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A24/raw_illumina_reads/A24_S2_R1_001.fastq

#in this run trim the Diagenode adapter error rate 20%
trim_galore -a GGTAGC -q 20 -e 0.2 --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC2/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore2 --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A24/raw_illumina_reads/A24_S2_R1_001.fastq


#I23
cd /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A26/
gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A26/raw_illumina_reads/A26_S3_R1_001.fastq.gz

#in this run trim the Diagenode adapter
trim_galore -q 20  --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A26/raw_illumina_reads/A26_S3_R1_001.fastq


#in this run trim the Diagenode adapter error rate 20%
trim_galore -a ATGAGC -q 20 -e 0.2 --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC2/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore2 --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A26/raw_illumina_reads/A26_S3_R1_001.fastq


#I24
cd /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A28/
gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A28/raw_illumina_reads/A28_S4_R1_001.fastq.gz

#in this run trim the Diagenode adapter
trim_galore -q 20  --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A28/raw_illumina_reads/A28_S4_R1_001.fastq

#in this run trim the Diagenode adapter error rate 20%
trim_galore -a CAAAAG -q 20 -e 0.2 --phred33 --gzip  --fastqc_args "--outdir /mnt/data1/Joana/MATRICS/RRBS_test/FASTQC2/" --output_dir /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore2 --rrbs  /mnt/data1/Joana/MATRICS/RRBS_test/Sample_A28/raw_illumina_reads/A28_S4_R1_001.fastq


#according to http://www.epigenesys.eu/images/stories/protocols/pdf/20120720103700_p57.pdf 'a deduplication step of the alignments is not desirable as it couuld remove a very large fraction of the aligned data altogether'


#check trimmed data on FASTQC


#run Bismark: NOTE assumes directional.

cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark
gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A17_S1_R1_001_trimmed.fq.gz 
 nohup bismark -N 1 --un --ambiguous --gzip -q --se /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A17_S1_R1_001_trimmed.fq --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1/ --genome_folder /mnt/data1/reference_files/mouse_GRCm38.p5/appended  > nohup_A17.out 2>&1&
 #29765

 
 #to run in the backgroup
 gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A24_S2_R1_001_trimmed.fq.gz 
 cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark
 nohup bismark -N 1 --un --ambiguous --gzip -q --se /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A24_S2_R1_001_trimmed.fq --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1/ --genome_folder /mnt/data1/reference_files/mouse_GRCm38.p5/appended > nohup_A24.out 2>&1&
 #30174
 

 
  #to run in the backgroup
 gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A26_S3_R1_001_trimmed.fq.gz 
 cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark
 nohup bismark -N 1 --un --ambiguous --gzip -q --se /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A26_S3_R1_001_trimmed.fq --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1/ --genome_folder /mnt/data1/reference_files/mouse_GRCm38.p5/appended > nohup_A26.out 2>&1&
 #31912
 
 
 
   #to run in the backgroup
 gzip -d  /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A28_S4_R1_001_trimmed.fq.gz 
 cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark
 nohup bismark -N 1 --un --ambiguous --gzip -q --se /mnt/data1/Joana/MATRICS/RRBS_test/trim_galore/A28_S4_R1_001_trimmed.fq --path_to_bowtie /mnt/data1/programs/bowtie2-2.3.1/ --genome_folder /mnt/data1/reference_files/mouse_GRCm38.p5/appended > nohup_A28.out 2>&1&
 #32016
 
 
 
 

 cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark/methylation_extractor
  nohup bismark_methylation_extractor  -s --gzip --comprehensive --multicore 4 --bedGraph /mnt/data1/Joana/MATRICS/RRBS_test/bismark/A17_S1_R1_001_trimmed_bismark_bt2.bam  > nohup_A17.out 2>&1&
49553
 
  cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark/methylation_extractor
  nohup bismark_methylation_extractor  -s --gzip --comprehensive --multicore 4 --bedGraph /mnt/data1/Joana/MATRICS/RRBS_test/bismark/A24_S2_R1_001_trimmed_bismark_bt2.bam  > nohup_A24.out 2>&1&
49578
 
  
  cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark/methylation_extractor
  nohup bismark_methylation_extractor  -s --gzip --comprehensive --multicore 2 --bedGraph /mnt/data1/Joana/MATRICS/RRBS_test/bismark/A26_S3_R1_001_trimmed_bismark_bt2.bam  > nohup_A26.out 2>&1&
49935
  
    
  cd /mnt/data1/Joana/MATRICS/RRBS_test/bismark/methylation_extractor
  nohup bismark_methylation_extractor  -s --gzip --comprehensive --multicore 2 --bedGraph /mnt/data1/Joana/MATRICS/RRBS_test/bismark/A28_S4_R1_001_trimmed_bismark_bt2.bam  > nohup_A28.out 2>&1&
49603
  
  
  