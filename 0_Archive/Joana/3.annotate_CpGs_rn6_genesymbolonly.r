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

read.csv("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/biomart_rat_270918.csv")->ensembl

gene.obj<-readTranscriptFeatures("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/ENS_Gene2_rn6_Annotation.bed",unique.prom=FALSE)

cpg.obj<-readFeatureFlank("/gpfs/ts0/projects/Research_Project-129726/MATRICS/general_files/UCSC_Gene2_rn6_GpGIslands_Annotation.bed", feature.flank.name=c("CpGi", "shores"))


RatExons<-as.data.frame(gene.obj$exons)
RatIntrons<-as.data.frame(gene.obj$introns)
RatPromoters<-as.data.frame(gene.obj$promoters)
RatTSSes<-as.data.frame(gene.obj$TSSes)

#keep only chromosomes 1 to 19, X and Y and MT
RatExons[which(as.character(RatExons$seqnames)%in%chroms), ]->RatExons_chr
RatIntrons[which(as.character(RatIntrons$seqnames)%in%chroms), ]->RatIntrons_chr
RatPromoters[which(as.character(RatPromoters$seqnames)%in%chroms), ]->RatPromoters_chr
RatTSSes[which(as.character(RatTSSes$seqnames)%in%chroms), ]->RatTSSes_chr


#add column with gene symbol
RatExons_chr$transcript<-gsub("\\..*","", RatExons_chr$name)
RatExons_chr$gene_symbol<-NA
RatExons_chr$gene_symbol<-ensembl$Gene.name[match(RatExons_chr$transcript, ensembl$Transcript.stable.ID)]

RatIntrons_chr$transcript<-gsub("\\..*","", RatIntrons_chr$name)
RatIntrons_chr$gene_symbol<-NA
RatIntrons_chr$gene_symbol<-ensembl$Gene.name[match(RatIntrons_chr$transcript, ensembl$Transcript.stable.ID)]

RatPromoters_chr$transcript<-gsub("\\..*","", RatPromoters_chr$name)
RatPromoters_chr$gene_symbol<-NA
RatPromoters_chr$gene_symbol<-ensembl$Gene.name[match(RatPromoters_chr$transcript, ensembl$Transcript.stable.ID)]

RatTSSes_chr$transcript<-gsub("\\..*","", RatTSSes_chr$name)
RatTSSes_chr$gene_symbol<-NA
RatTSSes_chr$gene_symbol<-ensembl$Gene.name[match(RatTSSes_chr$transcript, ensembl$Transcript.stable.ID)]


### create annotation frame for locations within 10kb of the start of each gene
RatGenes.ids<-unique(RatExons_chr$name)
RatGenes<-matrix(data = NA, nrow = length(RatGenes.ids), ncol = 4) 
for(i in 1:length(RatGenes.ids)){

	sub.exons<-RatExons_chr[which(RatExons_chr$name == RatGenes.ids[i]),]
	RatGenes[i,1]<-as.character(sub.exons$seqnames[1])
	RatGenes[i,2]<-min(sub.exons$start)
	RatGenes[i,3]<-max(sub.exons$end)
	RatGenes[i,4]<-RatGenes.ids[i]
}

## Add and take away 10kb from start and stop

RatGenes[,2]<-as.numeric(RatGenes[,2]) - 10000
RatGenes[,3]<-as.numeric(RatGenes[,3]) + 10000
RatGenes<-as.data.frame(RatGenes)
RatGenes[,2]<-as.numeric(as.character(RatGenes[,2]))
RatGenes[,3]<-as.numeric(as.character(RatGenes[,3]))
#RatGenes[,4]<-as.numeric(as.character(RatGenes[,4]))
colnames(RatGenes)<-c("seqnames", "start", "end", "name")

#add column with gene symbol
RatGenes$transcript<-gsub("\\..*","", RatGenes$name)
RatGenes$gene_symbol<-NA
RatGenes$gene_symbol<-NA
RatGenes$gene_symbol<-ensembl$Gene.name[match(RatGenes$transcript, ensembl$Transcript.stable.ID)]


RatCpGs<-as.data.frame(cpg.obj$CpGi)
RatCpGs[which(as.character(RatCpGs$seqnames)%in%chroms), ]->RatCpGs_chr

RatCpGShores<-as.data.frame(cpg.obj$shores)
RatCpGShores[which(as.character(RatCpGShores$seqnames)%in%chroms), ]->RatCpGShores_chr



GeneAnno<-vector(length = nrow(all_chr))
CpGAnno<-vector(length = nrow(all_chr))

library(doParallel)

cl<-makeCluster(8)
registerDoParallel(cl)

SearchRowGeneAnno<-function(i){
	
	output<-""
	position<- paste(all_chr$chrom[i], all_chr$chromStart[i], sep=":")
	
	### Check Exons
	tmp.sub<-subset(RatExons_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, " Exon ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 
	
	### Check Introns
	tmp.sub<-subset(RatIntrons_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "Intron ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 
	
	### Check Promoters
	tmp.sub<-subset(RatPromoters, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "Promoter ", tmp.sub$name[j], " - ", tmp.sub$gene_symbol[j], " | ")
		}
		
	} 

	### Check TSSs
	tmp.sub<-subset(RatTSSes, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(tmp.sub)){
			output<-paste0(output, "TSS ", tmp.sub$name[j], " | ")
		}
		
	} 

	if(output == ""){
		### Check within 10kb of gene
		tmp.sub<-subset(RatGenes, as.character(seqnames) == as.character(all_chr$chrom[i]))
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
	tmp.sub<-subset(RatCpGs_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
	tmp.sub<-subset(tmp.sub, start <= all_chr$chromStart[i])
	tmp.sub<-subset(tmp.sub, end >= all_chr$chromStart[i])
	
	if(nrow(tmp.sub) > 0){
		outputCpG<-paste(outputCpG, tmp.sub$name, "|")
		}
		
	### Check CpG shoress
	tmp.sub<-subset(RatCpGShores_chr, as.character(seqnames) == as.character(all_chr$chrom[i]))
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
