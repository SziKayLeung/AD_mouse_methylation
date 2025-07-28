setwd("/mnt/data1/Eilis/Projects/Rat_Neuron/") 
meth<-read.csv("CombinedMethylationValues_SummaryStats_minReadDepth10.csv", row.names = 1)


## For control pairs 1 and 2 rank abs diff multiple by sign and add ranks then take abs of rank

p1_rank<-rank(abs(meth$Diff_P1_NA))*sign(meth$Diff_P1_NA)
p2_rank<-rank(abs(meth$Diff_P2_NA))*sign(meth$Diff_P2_NA)
p3_rank<-rank(abs(meth$Diff_P3_FK_506))*sign(meth$Diff_P3_FK_506)
p4_rank<-rank(abs(meth$Diff_P4_nifedipine))*sign(meth$Diff_P4_nifedipine)
op_p3_rank<-max(p3_rank)-abs(p3_rank)+1
op_p4_rank<-max(p4_rank)-abs(p4_rank)+1


sum_rank_P1_P2<-abs(p1_rank+p2_rank)
sum_rank_P1_P2_P3<-abs(p1_rank+p2_rank+p3_rank)
sum_rank_P1_P2_P4<-abs(p1_rank+p2_rank+p4_rank)
sum_rank_P1_P2_oppP3<-abs(p1_rank+p2_rank)+op_p3_rank
sum_rank_P1_P2_oppP4<-abs(p1_rank+p2_rank)+op_p4_rank

meth<-cbind(meth, sum_rank_P1_P2,sum_rank_P1_P2_P3,sum_rank_P1_P2_P4, sum_rank_P1_P2_oppP3,sum_rank_P1_P2_oppP4)
meth<-meth[order(meth$sum_rank_P1_P2, decreasing = T), ]

###Filter so that min diff in both control pairs is as least 20%
sub<-subset(meth, abs(Diff_P1_NA) > 20)
sub<-subset(sub, abs(Diff_P2_NA) > 20)

### filter so that in same direction
sub<-subset(sub, Diff_P1_NA * Diff_P2_NA > 0)



### 

RatExons<-as.data.frame(gene.obj$exons)
RatIntrons<-as.data.frame(gene.obj$introns)
RatPromoters<-as.data.frame(gene.obj$promoters)
RatTSSes<-as.data.frame(gene.obj$TSSes)

RatExons$name<-gsub("RGD:", "", RatExons$name)
RatIntrons$name<-gsub("RGD:", "", RatIntrons$name)
RatPromoters$name<-gsub("RGD:", "", RatPromoters$name)
RatTSSes$name<-gsub("RGD:", "", RatTSSes$name)

RatCpGs<-as.data.frame(cpg.obj$CpGi)

GeneAnno<-vector(length = nrow(sub))
CpGAnno<-vector(length = nrow(sub))
for(i in 1:nrow(sub)){
	
	output<-""
	
	### Check Exons
	tmp.sub<-subset(RatExons, as.character(seqnames) == as.character(sub$chr[i]))
	tmp.sub<-subset(tmp.sub, start <= sub$chrpos[i])
	tmp.sub<-subset(tmp.sub, end >= sub$chrpos[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(rn4_geneanno.sub)){
			output<-paste(output, "Exon ", tmp.sub$name[j], "|")
		}
		
	} 
	
	### Check Introns
	tmp.sub<-subset(RatIntrons, as.character(seqnames) == as.character(sub$chr[i]))
	tmp.sub<-subset(tmp.sub, start <= sub$chrpos[i])
	tmp.sub<-subset(tmp.sub, end >= sub$chrpos[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(rn4_geneanno.sub)){
			output<-paste(output, "Intron ", tmp.sub$name[j], "|")
		}
		
	} 
	
	### Check Promoters
	tmp.sub<-subset(RatPromoters, as.character(seqnames) == as.character(sub$chr[i]))
	tmp.sub<-subset(tmp.sub, start <= sub$chrpos[i])
	tmp.sub<-subset(tmp.sub, end >= sub$chrpos[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(rn4_geneanno.sub)){
			output<-paste(output, "Promoter ", tmp.sub$name[j], "|")
		}
		
	} 

		### Check TSSs
	tmp.sub<-subset(RatTSSes, as.character(seqnames) == as.character(sub$chr[i]))
	tmp.sub<-subset(tmp.sub, start <= sub$chrpos[i])
	tmp.sub<-subset(tmp.sub, end >= sub$chrpos[i])
	if(nrow(tmp.sub) > 0){
		for(j in 1:nrow(rn4_geneanno.sub)){
			output<-paste(output, "TSS ", tmp.sub$name[j], "|")
		}
		
	} 
	
	GeneAnno[i]<-output
	
	### Check CpGs
	tmp.sub<-subset(RatCpGs, as.character(seqnames) == as.character(sub$chr[i]))
	tmp.sub<-subset(tmp.sub, start <= sub$chrpos[i])
	tmp.sub<-subset(tmp.sub, end >= sub$chrpos[i])
	
	if(nrow(tmp.sub) > 0){
		CpGAnno[i]<-tmp.sub$name
		}
}


sub<-cbind(sub, GeneAnno, CpGAnno)



write.csv(sub, "CpGsites_ConsistentDiffControlPairs_wGeneAnno.csv")

tab<-matrix(data = NA, ncol = 7, nrow = 1)
for(i in 1:7){
	sub1<-subset(sub, abs(Diff_P1_NA) > ((i*5)+15))
	tab[1,i]<-nrow(subset(sub1,  abs(Diff_P2_NA) > ((i*5)+15)))
	}

