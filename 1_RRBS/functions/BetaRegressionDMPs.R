## Isabel Castanho I.S.Castanho@exeter.ac.uk & Aisha Dahir A.N.Dahir@exeter.ac.uk
## Function for beta regression on RRBS complete dataset (for DMPs)
suppressMessages(library(betareg))
suppressMessages(library(lmtest)) # to use likelihood-ratio test (lrt)
suppressMessages(library(foreach))
suppressMessages(library(doParallel)) # set up parameters to run parallel

BetaRegressionRRBSinteraction <-  function(betas, pheno, formula, formulaNull, link = "probit", num.cores){
  # betas = dataframe containg betas
  # pheno = phenotypic file (coldata)
  # formula = full model (~Genotype + Age + Genotype*Age)
  # formulaNull = ~Genotype + Age
  # link = usually probit (same as Biseq)
  # no.cl = number of clusters to run parallel on
  
  ### set up parameters to run parallel
  cl<-makeCluster(num.cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(betareg))
  
  #Create function which performs analysis for each probe
  testCpG <- function(row, pheno, formula, formulaNull, link, min.meth, max.meth){
    
    library(betareg)
    library(lmtest)
    print(row)
    
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
    
    betas.j <- betas[j,]
    betas.j[betas.j == 0] <- min.meth
    betas.j[betas.j == 1] <- max.meth
    data <- cbind(t(betas.j), pheno)
    colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
    
    options(show.error.messages = FALSE)
    suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.null <- try(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])), data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel.null <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])), data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.lrt <- try(lrtest(stats.full, stats.null)))
    options(show.error.messages = TRUE)
    if((class(lmodel) == "try-error")){
      print("class(lmodel) is try-error")
    } else{
      if((class(stats.lrt)[1] == "try-error" | class(stats.lrt)[2] == "try-error" )){
        print("class(stats.lrt) is try-error")
      } else{
        if(!lmodel$converged){
          print("!lmodel$converged is TRUE")
        } else{
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          p.val.Genotype <- max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Genotype  #should not be zero
          baseline <- lmodel$coefficients$mean[1, 1]
          coef.Genotype <- lmodel$coefficients$mean[2, 1]
          se.Genotype <- lmodel$coefficients$mean[2, 2]
          meth.group1.WT <- inv.link(baseline) #meth.group1.WT
          meth.group2.TG <- inv.link(baseline + coef.Genotype) #meth.group2.TG
          meth.diff.Genotype <- as.numeric(meth.group1.WT) - as.numeric(meth.group2.TG)  #meth.diff.Genotype
          pseudo.R.sqrt <- lmodel$pseudo.r.squared #pseudo.R.sqrt
          estimate.Genotype <- coef.Genotype #estimate.Genotype
          std.error.Genotype <- se.Genotype #std.error.Genotype
          p.val.Age <- max(lmodel$coefficients$mean[3, 4], 1e-323) #p.val.Age      #should not be zero
          coef.Age <- lmodel$coefficients$mean[3, 1] 
          se.Age <- lmodel$coefficients$mean[3, 2] 
          estimate.Age <- coef.Age  #estimate.Age
          std.error.Age <- se.Age  #std.error.Age
          p.val.Interaction <- max(lmodel$coefficients$mean[4, 4], 1e-323) # p.val.Interaction # should not be zero
          coef.Interaction <- lmodel$coefficients$mean[4, 1] #estimate.Interaction
          se.Interaction <- lmodel$coefficients$mean[4, 2]  
          estimate.Interaction <- coef.Interaction #estimate.Interaction
          std.error.Interaction <- se.Interaction  #std.error.Interaction
          p.val.modelLRT <- max(stats.lrt$`Pr(>Chisq)`[2], 1e-323) #p.val.modelLRT
          Position <- rownames(betas.j)
          
          return(c(p.val.Genotype, meth.group1.WT, meth.group2.TG, meth.diff.Genotype,
             pseudo.R.sqrt, estimate.Genotype, std.error.Genotype,
             p.val.Age, estimate.Age, std.error.Age,
             p.val.Interaction, estimate.Interaction, std.error.Interaction, p.val.modelLRT,
             Position))
        }
      }
    }
    
    # return(c())
    
  }
  
  ###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
  min.meth <- min(betas[betas > 0], na.rm=TRUE)
  max.meth <- max(betas[betas < 1], na.rm=TRUE)
  
  res<-foreach(j=1:nrow(betas), .combine=rbind, .packages = c("betareg", "lmtest")) %dopar% {
    testCpG(betas[j,], pheno, formula = formula, formulaNull = formulaNull, link = link,
            max.meth = max.meth, min.meth = min.meth) 
  }
  
  stopCluster(cl)
  
  colnames(res) <- c("p.val.Genotype", "meth.group1.WT", "meth.group2.TG", "meth.diff.Genotype",
                     "pseudo.R.sqrt", "estimate.Genotype", "std.error.Genotype",
                     "p.val.Age", "estimate.Age", "std.error.Age",
                     "p.val.Interaction", "estimate.Interaction", "std.error.Interaction", "p.val.modelLRT",
                     "Position")
  
  rownames(res) <- rownames(betas)
  
  res<-as.data.frame(res)
  head(res)
  
  return(res)
}
