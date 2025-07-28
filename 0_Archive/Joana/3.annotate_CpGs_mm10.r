### Annotate all CpG sites found with >=1 read in all samples:

args <- commandArgs(TRUE) # option to read arguments from bash script
args[1]->cohort
args[2]->brain_region
args[3]->output_path
args[4]->genome

all<-read.table(paste0(output_path, "CpG_regions_txt_", cohort, brain_region, ".txt"), header=TRUE, stringsAsFactors=FALSE)

#rename mitochondrial chromosome
all$chrom[which(all$chrom=="chrMT")]<-"chrM"

#keep only chromosomes 1 to 19, X and Y and MT
paste0("chr", c(seq(1:19), "X", "Y", "M"))->chroms

all[which(all$chrom%in%chroms), ]->all_chr


### read in annotation files
library(methylKit) # works with module load R/3.3.3-intel-2017a-X11-20170314
library(GenomicRanges)
library(genomation)
#library(biomaRt)

#ensembl = useMart("ensembl")

#ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="www.ensembl.org")

#useDataset("mmusculus_gene_ensembl",mart=ensembl)

#check biomaRt has the latest version of the genome by doing:
#listDatasets(ensembl)[which(listDatasets(ensembl)[,1]=="mmusculus_gene_ensembl"),]

read.csv("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/biomart_mouse_270918.csv")->ensembl

gene.obj<-readTranscriptFeatures("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/ENS_Gene2_mm10_Annotation.bed",unique.prom=FALSE)

cpg.obj<-readFeatureFlank("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/UCSC_Gene2_mm10_GpGIslands_Annotation.bed", feature.flank.name=c("CpGi", "shores"))


MouseExons<-as.data.frame(gene.obj$exons)
MouseIntrons<-as.data.frame(gene.obj$introns)
MousePromoters<-as.data.frame(gene.obj$promoters)
MouseTSSes<-as.data.frame(gene.obj$TSSes)

#keep only chromosomes 1 to 19, X and Y and MT
MouseExons[which(as.character(MouseExons$seqnames)%in%chroms), ]->MouseExons_chr
MouseIntrons[which(as.character(MouseIntrons$seqnames)%in%chroms), ]->MouseIntrons_chr
MousePromoters[which(as.character(MousePromoters$seqnames)%in%chroms), ]->MousePromoters_chr
MouseTSSes[which(as.character(MouseTSSes$seqnames)%in%chroms), ]->MouseTSSes_chr


#add column with gene symbol
MouseExons_chr$transcript<-gsub("\\..*","", MouseExons_chr$name)
MouseExons_chr$gene_symbol<-NA
MouseExons_chr$gene_symbol<-ensembl$Gene.name[match(MouseExons_chr$transcript, ensembl$Transcript.stable.ID)]

MouseIntrons_chr$transcript<-gsub("\\..*","", MouseIntrons_chr$name)
MouseIntrons_chr$gene_symbol<-NA
MouseIntrons_chr$gene_symbol<-ensembl$Gene.name[match(MouseIntrons_chr$transcript, ensembl$Transcript.stable.ID)]

MousePromoters_chr$transcript<-gsub("\\..*","", MousePromoters_chr$name)
MousePromoters_chr$gene_symbol<-NA
MousePromoters_chr$gene_symbol<-ensembl$Gene.name[match(MousePromoters_chr$transcript, ensembl$Transcript.stable.ID)]

MouseTSSes_chr$transcript<-gsub("\\..*","", MouseTSSes_chr$name)
MouseTSSes_chr$gene_symbol<-NA
MouseTSSes_chr$gene_symbol<-ensembl$Gene.name[match(MouseTSSes_chr$transcript, ensembl$Transcript.stable.ID)]


### create annotation frame for locations within 10kb of the start of each gene
MouseGenes.ids<-unique(MouseExons_chr$name)
MouseGenes<-matrix(data = NA, nrow = length(MouseGenes.ids), ncol = 4) 
for(i in 1:length(MouseGenes.ids)){

	sub.exons<-MouseExons_chr[which(MouseExons_chr$name == MouseGenes.ids[i]),]
	MouseGenes[i,1]<-as.character(sub.exons$seqnames[1])
	MouseGenes[i,2]<-min(sub.exons$start)
	MouseGenes[i,3]<-max(sub.exons$end)
	MouseGenes[i,4]<-MouseGenes.ids[i]
}

## Add and take away 10kb from start and stop

MouseGenes[,2]<-as.numeric(MouseGenes[,2]) - 10000
MouseGenes[,3]<-as.numeric(MouseGenes[,3]) + 10000
MouseGenes<-as.data.frame(MouseGenes)
MouseGenes[,2]<-as.numeric(as.character(MouseGenes[,2]))
MouseGenes[,3]<-as.numeric(as.character(MouseGenes[,3]))
#MouseGenes[,4]<-as.numeric(as.character(MouseGenes[,4]))
colnames(MouseGenes)<-c("seqnames", "start", "end", "name")

#add column with gene symbol
MouseGenes$transcript<-gsub("\\..*","", MouseGenes$name)
MouseGenes$gene_symbol<-NA
MouseGenes$gene_symbol<-NA
MouseGenes$gene_symbol<-ensembl$Gene.name[match(MouseGenes$transcript, ensembl$Transcript.stable.ID)]


MouseCpGs<-as.data.frame(cpg.obj$CpGi)
MouseCpGs[which(as.character(MouseCpGs$seqnames)%in%chroms), ]->MouseCpGs_chr

MouseCpGShores<-as.data.frame(cpg.obj$shores)
MouseCpGShores[which(as.character(MouseCpGShores$seqnames)%in%chroms), ]->MouseCpGShores_chr



GeneAnno<-vector(length = nrow(all_chr))
CpGAnno<-vector(length = nrow(all_chr))

library(doParallel)

cl<-makeCluster(8)
registerDoParallel(cl)

SearchRowGeneAnno<-function(i){
	
	output<-""
	position<- paste(all_chr$chrom[i], all_chr$chromStart[i], sep=":")
	
	### Check Exons
	tmp.sub<-subset(MouseExons_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, " Exon ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 
	
	### Check Introns
	tmp.sub<-subset(MouseIntrons_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "Intron ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 
	
	### Check Promoters
	tmp.sub<-subset(MousePromoters, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "Promoter ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 

	### Check TSSs
	tmp.sub<-subset(MouseTSSes, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "TSS ", tmp.sub$name[j], " | ")
		}
		
	} 

	if(output == ""){
		### Check within 10kb of gene
		tmp.sub<-subset(MouseGenes, as.character(seqnames) == as.character(all_chr$chrom[i]))
		tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
		tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
		if(nrow(tmp.sub) > 0){
			dist<-vector(length = nrow(tmp.sub))
			for(j in 1:nrow(tmp.sub)){
				#output<-paste(output, "Within 10kb ", tmp.sub$name[j], "|")
				dist[j]<-min(c(abs(tmp.sub$start[j] + 10000 - all_chr$chromStart[i]), abs(tmp.sub$stop[j] - 10000 - all_chr$chromStart[i])))
			}
			output<-paste0(output, "Closest Gene Within 10 kb ", tmp.sub$name[which(dist == min(dist))], " - ", tmp.sub$gene_symbol[which(dist == min(dist))])
		
		}
	}
	
	return(c(position, output))
}

#GeneAnno<-foreach(i=1:nrow(all_chr), .combine = "c") %dopar% SearchRowGeneAnno(i)

GeneAnno<-matrix(NA, 1, 2)


for(i in 1:nrow(all_chr)){
SearchRowGeneAnno(i)->result
GeneAnno <- rbind(GeneAnno, result)
}

GeneAnno1 <- GeneAnno[-1,]

SearchRowCpGAnno<-function(i){	
	outputCpG<-""
	
	### Check CpGs
	tmp.sub<-subset(MouseCpGs_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	
	if(nrow(tmp.sub) > 0){
		outputCpG<-paste(outputCpG, tmp.sub$name, "|")
		}
		
	### Check CpG shoress
	tmp.sub<-subset(MouseCpGShores_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	
	if(nrow(tmp.sub) > 0){
		outputCpG<-paste(outputCpG, " shore |")
		}
		
	return(outputCpG)
}

CpGAnno<-foreach(i=1:nrow(all_chr), .combine = "c")%dopar% SearchRowCpGAnno(i)

all_chr1<-cbind(all_chr, GeneAnno1[,2], CpGAnno)

colnames(all_chr1) <- c(colnames(all_chr), "gene_annotation", "CpG_islands")

write.csv(all_chr1, paste0(output_path, cohort, "_", brain_region, "_Gene_Cpgislands_Anno.csv"))
