# Correlation matrix of RRBS samples against each other

# Dorothea Seiler Vellame 23-10-2018

# load the data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixFiltered.Rdata")
rTg<-RRBS
# keep methyl cols only
rTgm<-rTg[,2:63]

load("J20MatrixFiltered.Rdata")
J20<-RRBS
# keep methyl cols only
J20m<-J20[,2:64]
rm(RRBS)

# practice plot
plot(J20m[,1],J20m[,3],xlim=c(0,100),ylim=c(0,100))




df <- data.frame(J20m[,1],J20m[,3])

## Use densCols() output to get density at each point
x <- densCols(J20m[,1],J20m[,3], colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]

## Plot it, reordering rows so that densest points are plotted on top
plot(J20m[,1]~J20m[,3], data=df[order(df$dens),], pch=20, col=col, cex=2)



library(LSD)
## different plot test
heatscatter(J20m[,1],J20m[,2])
