# Function to plot methylation on section of genome with differential methylation returned from bumphunter

## INPUT: tab - output from bumphunter, tab is col 1 of tab$table, start and end are 2 and 3 (will be needed to plot)
##        chr and pos
##        mModel = J20 or rTg4510 - Tells which mouse model is used 
##        phenotype - contains phenotype info on the model
##        RRBS - with all NAs removed already (was done in bumphunter stages)
##        numberOfPlots = 1 or 2, whether you get 1 or 2 plots per graph
##        pdfFileName = dmrPlots.pdf

## OUTPUT: pdf containing plots of top dm regions either with 1 plot of all or 2 stacked plots


setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")

load("rTg4510BHRestult.Rdata")
load("posChrfile.Rdata")

#for testing 
r=1 # i refers to the row in tab
mModel="rTg4510"
numberOfPlots = 1
tabstore<-tab
tabsub=lapply(tab,head,n=100)
tab<-tabsub


bumphunterPlot(tab, chr, pos, RRBS, phenotype, mModel="rTg4510", numberOfPlots = 1)

bumphunterPlot = function (tab, chr, pos, RRBS, phenotype, mModel="rTg4510", numberOfPlots = 1, pdfFileName="dmrPlots.pdf"){
  for (r in 2:nrow(tab$table)){
    # isolate region to be plotted
    
    pos<-as.numeric(pos)
    print(r)
    print(tab$table[r,])
    k<-which(chr==tab$table[r,1])
    j<-which(pos>tab$table[r,2])
    l<-which(pos<tab$table[r,3])
    
    
    print(c(min(pos[j]),max(pos[j])))
    print(c(min(pos[l]),max(pos[l])))
    
    
    
    n<-intersect(intersect(k,j),l)
    
    print(n)
    
    if (length(n)==0){
      print(paste("gap size =",tab$table[r,3]-tab$table[r,2]))
      print(r)
    }
    
    # get methylation values in genomic region
    if (mModel == "J20"){
      RRBSplot<-cbind(as.numeric(pos[n]),RRBS[n,2:64])
    }
    
    if (mModel == "rTg4510"){
      RRBSplot<-cbind(as.numeric(pos[n]),RRBS[n,2:63])
    }
    
    
    # create numerical phenotype numbers 
    # 1=TG2 2=TG4 3=TG6 4=TG8 5=WT2 6=WT4 7=WT6 8=WT8  
    # months will depend on the model
    # NOTE TG before WT
    pheno<-as.character(phenotype[,6]) 
    j=1
    for (i in 1:length(levels(phenotype$Group_ID))){
      pheno[which(pheno==as.character(levels(phenotype$Group_ID)[i]))]<-j
      j=j+1
    }
    
    pheno<-as.numeric(pheno)
    
    library(RColorBrewer)
    colWT<-brewer.pal(4,"Greys")
    
    if (mModel == "J20"){
      colJ20<-brewer.pal(6,"OrRd")
      col<-c(colJ20[3:6],colWT)
    }
    if(mModel == "rTg4510"){
      colrTg<-brewer.pal(6,"YlGnBu")
      col<-c(colrTg[3:6],colWT)
    }
    
    pdf(file=pdfFileName)
    if (numberOfPlots==1){
      par(mfrow=c(1,1))
      plot(1, type="n", xlim=c(min(RRBSplot[,1]),max(RRBSplot[,1])),ylim=c(0,100), xlab=paste("Position on Chromosome",chr[n[1]]), ylab="% methylation")
      for (i in 2:ncol(RRBSplot)){
        points(RRBSplot[,1],RRBSplot[,i],col="black", bg=col[pheno[i-1]],pch=21,cex=0.9)
      }
      if (mModel=="rTg4510"){
        legend("bottomright",col="black", pt.bg=col,pch=21,cex=0.9,legend=c("TG2","TG4", "TG6", "TG8", "WT2", "WT4", "WT6", "WT8"),
               bty="n",horiz=TRUE,xpd=TRUE,inset=c(0,1))
      }
      if (mModel=="J20"){
        legend("bottomright",col="black", pt.bg=col,pch=21,cex=0.9,legend=c("TG6","TG8", "TG10", "TG12", "WT6", "WT8", "WT10", "WT12"),
               bty="n",horiz=TRUE,xpd=TRUE,inset=c(0,1))
      }
    }
    
    
    
    if (numberOfPlots==2){
      par(mfrow=c(2,1))
      plot(1, type="n", xlim=c(min(RRBSplot[,1]),max(RRBSplot[,1])),ylim=c(0,100), xlab=paste("Position on Chromosome",chr[n[1]]), ylab="% methylation")
      RRBSWT<-cbind(RRBSplot[,1],RRBSplot[,(which(phenotype$Genotype=="WT")+1)])
      RRBSTG<-cbind(RRBSplot[,1],RRBSplot[,(which(phenotype$Genotype=="TG")+1)])
      WTpheno<-pheno[(which(phenotype$Genotype=="WT"))]
      TGpheno<-pheno[(which(phenotype$Genotype=="TG"))]
      
      for (i in 2:ncol(RRBSWT)){
        points(RRBSplot[,1],RRBSWT[,i],col="black", bg=col[WTpheno[i-1]],pch=21)
      }
      if (mModel=="rTg4510"){
        legend("bottomright",col="black", pt.bg=col,pch=21,cex=0.9,legend=c("TG2","TG4", "TG6", "TG8", "WT2", "WT4", "WT6", "WT8"),
               bty="n",horiz=TRUE,xpd=TRUE,inset=c(0,1))
      }
      if (mModel=="J20"){
        legend("bottomright",col="black", pt.bg=col,pch=21,cex=0.9,legend=c("TG6","TG8", "TG10", "TG12", "WT6", "WT8", "WT10", "WT12"),
               bty="n",horiz=TRUE,xpd=TRUE,inset=c(0,1))
      }
      plot(1, type="n", xlim=c(min(RRBSplot[,1]),max(RRBSplot[,1])),ylim=c(0,100), xlab=paste("Position on Chromosome",chr[n[1]]), ylab="% methylation")
      for (i in 2:ncol(RRBSTG)){
        points(RRBSplot[,1],RRBSTG[,i],col="black", bg=col[TGpheno[i-1]],pch=21)
      }
    }
    dev.off()
  }
}
