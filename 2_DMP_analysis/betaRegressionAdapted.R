# Adapted Beta Regression function from Biseq (https://rdrr.io/bioc/BiSeq/src/R/betaRegression.R)

# Isabel Castanho (I.S.Castanho@exeter.ac.uk) and Aisha Dahir (A.N.Dahir@exeter.ac.uk)

library(BiSeq)
library(betareg)
library(lmtest) # to use likelihood-ratio test (lrt)

betaRegressionAdapted <- function(formula, formulaNull, link = "probit", object, mc.cores, ...){
  
  strand(object) <- "*"
  object <- sort(object)
  
  if(link == "loglog"){
    inv.link <- function(x){
      exp(-exp(-x))
    }
  }
  if(link == "logit"){
    inv.link <- function(x){
      1 / ( 1 + exp(-x))
    }
  }
  if(link == "probit"){
    inv.link <- function(x){
      pnorm(x)
    }
  }
  if(link == "cloglog"){
    inv.link <- function(x){
      1 - exp(-exp(x))
    }
  }
  if(link == "log"){
    inv.link <- function(x){
      exp(x)
    }
  }
  
  min.meth <- min(methLevel(object)[methLevel(object) > 0], na.rm=TRUE)
  max.meth <- max(methLevel(object)[methLevel(object) < 1], na.rm=TRUE)
  
  object.split <- split(object, f = rep(1:mc.cores, length.out = nrow(object)))
  
  beta.regr <- function(object.part, formula, formulaNull, link, ...){
    chr <- as.character(seqnames(rowRanges(object.part)))
    pred.meth.part <- methLevel(object.part)
    
    p.val.Genotype <- rep(NA,nrow(object.part))
    meth.group1.WT <- rep(NA,nrow(object.part))
    meth.group2.TG <- rep(NA,nrow(object.part))
    meth.diff.Genotype <- rep(NA,nrow(object.part))
    direction.Genotype <- rep(NA,nrow(object.part))
    pseudo.R.sqrt <- rep(NA,nrow(object.part))
    estimate.Genotype <- rep(NA,nrow(object.part))
    std.error.Genotype <- rep(NA,nrow(object.part))
    p.val.Age <- rep(NA,nrow(object.part))
    direction.Age <- rep(NA,nrow(object.part))
    estimate.Age <- rep(NA,nrow(object.part))
    std.error.Age <- rep(NA,nrow(object.part))
    p.val.Interaction <- rep(NA,nrow(object.part))
    direction.Interaction <- rep(NA,nrow(object.part))
    estimate.Interaction <- rep(NA,nrow(object.part))
    std.error.Interaction <- rep(NA,nrow(object.part))
    p.val.modelLRT <- rep(NA,nrow(object.part))
    
    
    for(j in 1:nrow(pred.meth.part)){
      pred.meth <- pred.meth.part[j,]
      pred.meth[pred.meth == 0] <- min.meth
      pred.meth[pred.meth == 1] <- max.meth
      data <- cbind(pred.meth = pred.meth,
                    as.data.frame(colData(object)))
      options(show.error.messages = FALSE)
      suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("pred.meth ~ ", formula[2])), data=data, link = link, ...), silent=TRUE))
      suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("pred.meth ~ ", formula[2])),data=data, link = link, ...)), silent=TRUE))
      suppressWarnings(stats.null <- try(betareg(formula = as.formula(paste("pred.meth ~ ", formulaNull[2])),data=data, link = link, ...), silent=TRUE))
      suppressWarnings(lmodel.null <- try(summary(betareg(formula = as.formula(paste("pred.meth ~ ", formulaNull[2])),data=data, link = link, ...)), silent=TRUE))
      suppressWarnings(stats.lrt <- try(lrtest(stats.full, stats.null)))
      options(show.error.messages = TRUE)
      if((class(lmodel) == "try-error")){
        print("class(lmodel) is try-error")
        p.val.Genotype[j] <- NA
        meth.diff.Genotype[j] <- NA
        p.val.Age[j] <- NA
        p.val.Interaction[j] <- NA
      } else{
        if((class(stats.lrt) == "try-error")){
          print("class(stats.lrt) is try-error")
          p.val.modelLRT[j] <- NA
        } else{
          if(!lmodel$converged){
            print("!lmodel$converged is TRUE")
            p.val.Genotype[j] <- NA
            meth.diff.Genotype[j] <- NA
            p.val.Age[j] <- NA
            p.val.Interaction[j] <- NA
          } else{
            # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
            # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
            p.val.Genotype[j] <- max(lmodel$coefficients$mean[2, 4], 1e-323) # should not be zero
            baseline <- lmodel$coefficients$mean[1, 1]
            coef.Genotype <- lmodel$coefficients$mean[2, 1]
            se.Genotype <- lmodel$coefficients$mean[2, 2]
            meth.group1.WT[j] <- inv.link(baseline)
            meth.group2.TG[j] <- inv.link(baseline + coef.Genotype)
            meth.diff.Genotype[j] <-  meth.group1.WT[j] - meth.group2.TG[j]
            pseudo.R.sqrt[j] <- lmodel$pseudo.r.squared
            estimate.Genotype[j] <- coef.Genotype
            std.error.Genotype[j] <- se.Genotype
            p.val.Age[j] <- max(lmodel$coefficients$mean[3, 4], 1e-323) # should not be zero
            coef.Age <- lmodel$coefficients$mean[3, 1]
            se.Age <- lmodel$coefficients$mean[3, 2]
            estimate.Age[j] <- coef.Age
            std.error.Age[j] <- se.Age
            p.val.Interaction[j] <- max(lmodel$coefficients$mean[4, 4], 1e-323) # should not be zero
            coef.Interaction <- lmodel$coefficients$mean[4, 1]
            se.Interaction <- lmodel$coefficients$mean[4, 2]
            estimate.Interaction[j] <- coef.Interaction
            std.error.Interaction[j] <- se.Interaction
            p.val.modelLRT[j] <- max(stats.lrt$`Pr(>Chisq)`[2], 1e-323)
          }
        }
      }
    }
    
    out <- data.frame(chr = chr,
                      pos = start(ranges(object.part)),
                      p.val.Genotype,
                      meth.group1.WT,
                      meth.group2.TG,
                      meth.diff.Genotype,
                      estimate.Genotype,
                      std.error.Genotype,
                      pseudo.R.sqrt,
                      p.val.Age,
                      estimate.Age,
                      std.error.Age,
                      p.val.Interaction,
                      estimate.Interaction,
                      std.error.Interaction,
                      p.val.modelLRT)
    print(head(out))
    return(out)
  }
  
  summary.df.l <- mclapply(object.split, beta.regr, formula=formula, formulaNull=formulaNull, link=link, mc.cores=mc.cores, ...)
  summary.df <- do.call(rbind, summary.df.l)
  # summary.df <- as.data.frame(summary.df) # couldnt use $ in matrix so converted to a dataframe
  pos.object <- paste(seqnames(object), start(object), sep="_")
  pos.summary <- paste(summary.df$chr, summary.df$pos, sep="_")
  ind <- match(pos.object, pos.summary)
  summary.df <- summary.df[ind,]
  summary.df$cluster.id <- mcols(rowRanges(object))$cluster.id
  return(summary.df)
}



#setMethod("betaRegressionAdapted",
#          signature=c(formula = "formula", formulaNull = "formula", link = "character", object="BSrel", mc.cores = "numeric"),
#          .betaRegressionAdapted)

#setMethod("betaRegressionAdapted",
#          signature=c(formula = "formula", formulaNull = "formula", link = "character",  object="BSrel", mc.cores = "missing"),
#          function(formula,formulaNull, link, object, ...) {
#            .betaRegressionAdapted(formula, formulaNull, link, object, mc.cores = 1, ...)
#          })