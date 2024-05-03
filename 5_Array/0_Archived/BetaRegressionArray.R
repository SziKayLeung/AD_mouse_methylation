## Aisha Dahir A.N.Dahir@exeter.ac.uk
library(betareg)
library(lmtest) # to use likelihood-ratio test (lrt)
## Function for beta regression on array dataset
BetaRegressionArray <-  function(betas, pheno, formula, formulaNull, link = "probit"){
  
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
  
  min.meth <- min(betas[betas > 0], na.rm=TRUE)
  max.meth <- max(betas[betas < 1], na.rm=TRUE)
  
  res <- matrix(NA, ncol = 15, nrow = nrow(betas))
  colnames(res) <- c("p.val.Genotype","meth.group1.WT",
                     "meth.group2.TG", "meth.diff.Genotype", 
                     "pseudo.R.sqrt","estimate.Genotype","std.error.Genotype", "p.val.Age", 
                     "estimate.Age", "std.error.Age", "p.val.Interaction",
                     "estimate.Interaction", "std.error.Interaction",
                     "p.val.modelLRT", "cpg")
  
  for(j in 1:nrow(betas)){
    betas.j <- betas[j,]
    betas.j[betas.j == 0] <- min.meth
    betas.j[betas.j == 1] <- max.meth
    data <- cbind(t(betas.j), pheno)
    colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
    
    
    options(show.error.messages = FALSE)
    suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.null <- try(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel.null <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.lrt <- try(lrtest(stats.full, stats.null)))
    options(show.error.messages = TRUE)
    if((class(lmodel) == "try-error")){
      print("class(lmodel) is try-error")
      res[j,1:14] <- rep(NA, 14)
    } else{
      if((class(stats.lrt)[1] == "try-error" | class(stats.lrt)[2] == "try-error" )){
        print("class(stats.lrt) is try-error")
        res[j,1:14] <- rep(NA, 14)
      } else{
        if(!lmodel$converged){
          print("!lmodel$converged is TRUE")
          res[j,1:14] <- rep(NA, 14)
        } else{
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          res[j,1] <- max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Genotype  #should not be zero
          baseline <- lmodel$coefficients$mean[1, 1]
          coef.Genotype <- lmodel$coefficients$mean[2, 1]
          se.Genotype <- lmodel$coefficients$mean[2, 2]
          res[j,2] <- inv.link(baseline) #meth.group1.WT
          res[j,3] <- inv.link(baseline + coef.Genotype) #meth.group1.TG
          res[j,4] <- as.numeric(res[j,2]) - as.numeric(res[j,3])  #meth.diff.Genotype
          res[j,5] <- lmodel$pseudo.r.squared #pseudo.R.sqrt
          res[j,6] <- coef.Genotype #estimate.Genotype
          res[j,7] <- se.Genotype #std.error.Genotype
          res[j,8] <- max(lmodel$coefficients$mean[3, 4], 1e-323) #p.val.Age      #should not be zero
          coef.Age <- lmodel$coefficients$mean[3, 1] 
          se.Age <- lmodel$coefficients$mean[3, 2] 
          res[j,9] <- coef.Age  #estimate.Age
          res[j,10] <- se.Age  #std.error.Age
          res[j,11] <- max(lmodel$coefficients$mean[7, 4], 1e-323) # p.val.Interaction # should not be zero
          coef.Interaction <- lmodel$coefficients$mean[7, 1] #estimate.Interaction
          se.Interaction <- lmodel$coefficients$mean[7, 2]  
          res[j,12] <- coef.Interaction #estimate.Interaction
          res[j,13] <- se.Interaction  #std.error.Interaction
          res[j,14] <- max(stats.lrt$`Pr(>Chisq)`[2], 1e-323) #p.val.modelLRT
          
          res[,15] <- rownames(betas)
          
        }
      }
    }
  }
  return(res)
}


BetaRegressionFullModel <-  function(betas, pheno, formula, formulaNull, link = "probit"){
  
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
  
  min.meth <- min(betas[betas > 0], na.rm=TRUE)
  max.meth <- max(betas[betas < 1], na.rm=TRUE)
  
  res <- matrix(NA, ncol = 14, nrow = nrow(betas))
  colnames(res) <- c("p.val.Genotype","meth.group1.WT",
                     "meth.group2.TG", "meth.diff.Genotype", 
                     "pseudo.R.sqrt","estimate.Genotype","std.error.Genotype", "p.val.Age", 
                     "estimate.Age", "std.error.Age", "p.val.Pathology",
                     "estimate.Pathology", "std.error.Pathology", "cpg")
  
  for(j in 1:nrow(betas)){
    betas.j <- betas[j,]
    betas.j[betas.j == 0] <- min.meth
    betas.j[betas.j == 1] <- max.meth
    data <- cbind(t(betas.j), pheno)
    colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
    options(show.error.messages = FALSE)
    suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),data=data, link = link)), silent=TRUE))
    options(show.error.messages = TRUE)
    if((class(lmodel) == "try-error")){
      print("class(lmodel) is try-error")
      res[j,1:14] <- rep(NA, 14)
    } else{
      if(!lmodel$converged){
        print("!lmodel$converged is TRUE")
        res[j,1:14] <- rep(NA, 14)
      } else{
        # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
        # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
        res[j,1] <- max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Genotype  #should not be zero
        baseline <- lmodel$coefficients$mean[1, 1]
        coef.Genotype <- lmodel$coefficients$mean[2, 1]
        se.Genotype <- lmodel$coefficients$mean[2, 2]
        res[j,2] <- inv.link(baseline) #meth.group1.WT
        res[j,3] <- inv.link(baseline + coef.Genotype) #meth.group1.TG
        res[j,4] <- as.numeric(res[j,2]) - as.numeric(res[j,3])  #meth.diff.Genotype
        res[j,5] <- lmodel$pseudo.r.squared #pseudo.R.sqrt
        res[j,6] <- coef.Genotype #estimate.Genotype
        res[j,7] <- se.Genotype #std.error.Genotype
        res[j,8] <- max(lmodel$coefficients$mean[3, 4], 1e-323) #p.val.Age      #should not be zero
        coef.Age <- lmodel$coefficients$mean[3, 1] 
        se.Age <- lmodel$coefficients$mean[3, 2] 
        res[j,9] <- coef.Age  #estimate.Age
        res[j,10] <- se.Age  #std.error.Age
        res[j,11] <- max(lmodel$coefficients$mean[7, 4], 1e-323) # p.val.Pathology # should not be zero
        coef.Pathology <- lmodel$coefficients$mean[7, 1] #estimate.Pathology
        se.Pathology <- lmodel$coefficients$mean[7, 2]  
        res[j,12] <- coef.Pathology #estimate.Pathology
        res[j,13] <- se.Pathology  #std.error.Pathology
        
        res[,14] <- rownames(betas)
        
      }
    }
  }
  return(res)
}



BetaRegressionArrayPathology <-  function(betas, pheno, formula, link = "probit"){
  
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
  
  min.meth <- min(betas[betas > 0], na.rm=TRUE)
  max.meth <- max(betas[betas < 1], na.rm=TRUE)
  
  res <- matrix(NA, ncol = 8, nrow = nrow(betas))
  colnames(res) <- c("p.val.Pathology","meth.group1",
                     "meth.group2", "meth.diff", 
                     "pseudo.R.sqrt","estimate.Pathology",
                     "std.error.Pathology", "cpg")
  
  for(j in 1:nrow(betas)){
    betas.j <- betas[j,]
    betas.j[betas.j == 0] <- min.meth
    betas.j[betas.j == 1] <- max.meth
    data <- cbind(t(betas.j), pheno)
    colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
    
    
    options(show.error.messages = FALSE)
    suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),data=data, link = link)), silent=TRUE))
    options(show.error.messages = TRUE)
    if((class(lmodel) == "try-error")){
      print("class(lmodel) is try-error")
      res[j,1:14] <- rep(NA, 14)
    } else{
        if(!lmodel$converged){
          print("!lmodel$converged is TRUE")
          res[j,1:14] <- rep(NA, 14)
        } else{
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          res[j,1] <- max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Genotype  #should not be zero
          baseline <- lmodel$coefficients$mean[1, 1]
          coef.Genotype <- lmodel$coefficients$mean[2, 1]
          se.Genotype <- lmodel$coefficients$mean[2, 2]
          res[j,2] <- inv.link(baseline) #meth.group1.WT
          res[j,3] <- inv.link(baseline + coef.Genotype) #meth.group1.TG
          res[j,4] <- as.numeric(res[j,2]) - as.numeric(res[j,3])  #meth.diff.Genotype
          res[j,5] <- lmodel$pseudo.r.squared #pseudo.R.sqrt
          res[j,6] <- coef.Genotype #estimate.Genotype
          res[j,7] <- se.Genotype #std.error.Genotype
          res[,8] <- rownames(betas)
          
        }
      }
    }
  return(res)
}

BetaRegGenotypeNoChip <-  function(betas, pheno, formula, formulaNull, link = "probit"){
  
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
  
  min.meth <- min(betas[betas > 0], na.rm=TRUE)
  max.meth <- max(betas[betas < 1], na.rm=TRUE)
  
  res <- matrix(NA, ncol = 15, nrow = nrow(betas))
  colnames(res) <- c("p.val.Genotype","meth.group1.WT",
                     "meth.group2.TG", "meth.diff.Genotype", 
                     "pseudo.R.sqrt","estimate.Genotype","std.error.Genotype", "p.val.Age", 
                     "estimate.Age", "std.error.Age", "p.val.Interaction",
                     "estimate.Interaction", "std.error.Interaction",
                     "p.val.modelLRT", "cpg")
  
  for(j in 1:nrow(betas)){
    betas.j <- betas[j,]
    betas.j[betas.j == 0] <- min.meth
    betas.j[betas.j == 1] <- max.meth
    data <- cbind(t(betas.j), pheno)
    colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
    
    
    options(show.error.messages = FALSE)
    suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.null <- try(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link), silent=TRUE))
    suppressWarnings(lmodel.null <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link)), silent=TRUE))
    suppressWarnings(stats.lrt <- try(lrtest(stats.full, stats.null)))
    options(show.error.messages = TRUE)
    if((class(lmodel) == "try-error")){
      print("class(lmodel) is try-error")
      res[j,1:14] <- rep(NA, 14)
    } else{
      if((class(stats.lrt)[1] == "try-error" | class(stats.lrt)[2] == "try-error" )){
        print("class(stats.lrt) is try-error")
        res[j,1:14] <- rep(NA, 14)
      } else{
        if(!lmodel$converged){
          print("!lmodel$converged is TRUE")
          res[j,1:14] <- rep(NA, 14)
        } else{
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          res[j,1] <- max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Genotype  #should not be zero
          baseline <- lmodel$coefficients$mean[1, 1]
          coef.Genotype <- lmodel$coefficients$mean[2, 1]
          se.Genotype <- lmodel$coefficients$mean[2, 2]
          res[j,2] <- inv.link(baseline) #meth.group1.WT
          res[j,3] <- inv.link(baseline + coef.Genotype) #meth.group1.TG
          res[j,4] <- as.numeric(res[j,2]) - as.numeric(res[j,3])  #meth.diff.Genotype
          res[j,5] <- lmodel$pseudo.r.squared #pseudo.R.sqrt
          res[j,6] <- coef.Genotype #estimate.Genotype
          res[j,7] <- se.Genotype #std.error.Genotype
          res[j,8] <- max(lmodel$coefficients$mean[3, 4], 1e-323) #p.val.Age      #should not be zero
          coef.Age <- lmodel$coefficients$mean[3, 1] 
          se.Age <- lmodel$coefficients$mean[3, 2] 
          res[j,9] <- coef.Age  #estimate.Age
          res[j,10] <- se.Age  #std.error.Age
          res[j,11] <- max(lmodel$coefficients$mean[4, 4], 1e-323) # p.val.Interaction # should not be zero
          coef.Interaction <- lmodel$coefficients$mean[4, 1] #estimate.Interaction
          se.Interaction <- lmodel$coefficients$mean[4, 2]  
          res[j,12] <- coef.Interaction #estimate.Interaction
          res[j,13] <- se.Interaction  #std.error.Interaction
          res[j,14] <- max(stats.lrt$`Pr(>Chisq)`[2], 1e-323) #p.val.modelLRT
          
          res[,15] <- rownames(betas)
          
        }
      }
    }
  }
  return(res)
}



#### Running in parallel


  
testCpG<-function(row, pheno, max.meth, min.meth, link, formula, formulaNull){
  
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
  
  betas.j <- row
  betas.j[betas.j == 0] <- min.meth
  betas.j[betas.j == 1] <- max.meth
  data <- cbind(t(betas.j), pheno)
  colnames(data)[1] <- "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
  
  options(show.error.messages = FALSE)
  suppressWarnings(stats.full <- try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])), data=data, link = link), silent=TRUE))
  suppressWarnings(lmodel <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),data=data, link = link)), silent=TRUE))
  suppressWarnings(stats.null <- try(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link), silent=TRUE))
  suppressWarnings(lmodel.null <- try(summary(betareg(formula = as.formula(paste("betas.j ~ ", formulaNull[2])),data=data, link = link)), silent=TRUE))
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
        meth.group2.TG <- inv.link(baseline + coef.Genotype) #meth.group1.TG
        meth.diff.Genotype <- as.numeric(res[j,2]) - as.numeric(res[j,3])  #meth.diff.Genotype
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
        
      }
    }
  }
  return(c(p.val.Genotype, meth.group1.WT, meth.group2.TG, meth.diff.Genotype, pseudo.R.sqrt, estimate.Genotype, std.error.Genotype,
           p.val.Age, estimate.Age, std.error.Age, p.val.Interaction, estimate.Interaction, std.error.Interaction, p.val.modelLRT ))
}


