# Create a function that makes box plots and same with lines instead per all sample covs passed in

# INPUT: RRBScov - the coverage columns of the RRBS data
#        box=TRUE - Plot boxplots?
#        lines=TRUE - Plot line box plot?   One of box or lines should be true or this function wont do anything

# Dorothea Seiler Vellame 26-10-2018

# get test input
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("J20MatrixUnfiltered.Rdata")
RRBScov<-RRBSmatrix[,65:127]


covBoxPlots<-function(RRBScov,box=TRUE,lines=TRUE){
  
  boxplot(RRBScov,log="y",ylab="log coverage")
  
  
  
}

boxplot(RRBScov[,1:3],log="y",ylab="log coverage")


# Example data
tmin <- as.Date("2000-01-01")
tmax <- as.Date("2001-01-01")
tlab <- seq(tmin, tmax, by="month")
lab <- format(tlab,format="%Y-%b")
set.seed(111)
x <- seq(tmin, tmax,100)
y <- cumsum(rnorm(100))

# Plot
png("plot_w_rotated_axis_labels.png", height=3, width=6, units="in", res=400)
op <- par(mar=c(6,4,1,1))
plot(x, y, t="l", xaxt="n", xlab="")
axis(1, at=tlab, labels=FALSE)
text(x=tlab, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels=lab, srt=90, adj=1, xpd=TRUE)
par(op)
dev.off()

ltab<-tlab[1:3]

boxplot(RRBScov[,1:3],log="y",ylab="log coverage", xaxt="n", xlab="")
axis(1, at=c(1,2,3), labels=FALSE)
text(x=c(1,2,3), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels=lab, srt=90, adj=1, xpd=TRUE)
