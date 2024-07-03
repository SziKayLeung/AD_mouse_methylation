## Isabel Castanho I.S.Castanho@exeter.ac.uk
## Function for beta regression on RRBS complete dataset (for DMPs)
library(betareg)
library(lmtest) # to use likelihood-ratio test (lrt)
library(foreach)
library(doParallel) # set up parameters to run parallel

BetaRegressionRRBSpathology <-
  function(betas, pheno, formula, link = "probit", num.cores) {
    # betas = dataframe containg betas
    # pheno = phenotypic file (coldata)
    # formula = ~ECX (Pathology)
    # link = usually probit (same as Biseq)
    # no.cl = number of clusters to run parallel on
    
    ### set up parameters to run parallel
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)
    clusterEvalQ(cl, library(betareg))
    
    #Create function which performs analysis for each probe
    testCpG <-
      function(row,
               pheno,
               formula,
               link,
               min.meth,
               max.meth) {
        # library(betareg)
        # library(lmtest)
        
        if (link == "loglog") {
          inv.link <- function(x) {
            exp(-exp(-x))
          }
        }
        if (link == "logit") {
          inv.link <- function(x) {
            1 / (1 + exp(-x))
          }
        }
        if (link == "probit") {
          inv.link <- function(x) {
            pnorm(x)
          }
        }
        if (link == "cloglog") {
          inv.link <- function(x) {
            1 - exp(-exp(x))
          }
        }
        if (link == "log") {
          inv.link <- function(x) {
            exp(x)
          }
        }
        
        betas.j <- betas[j, ]
        betas.j[betas.j == 0] <- min.meth
        betas.j[betas.j == 1] <- max.meth
        data <- cbind(t(betas.j), pheno)
        colnames(data)[1] <-
          "betas.j" #usually puts probe name but we need to keep the header consistent - so go with beta.j
        
        options(show.error.messages = FALSE)
        suppressWarnings(stats.full <-
                           try(betareg(formula = as.formula(paste("betas.j ~ ", formula[2])),
                                       data = data,
                                       link = link),
                               silent = TRUE)
        )
        suppressWarnings(lmodel <-
                           try(summary(
                             betareg(
                               formula = as.formula(paste("betas.j ~ ", formula[2])),
                               data  =  data,
                               link = link
                             )
                           )
                           , silent  =  TRUE)
        )
        options(show.error.messages = TRUE)
        if ((class(lmodel) == "try-error")) {
          print("class(lmodel) is try-error")
        } else {
          if (!lmodel$converged) {
            print("!lmodel$converged is TRUE")
          } else {
            # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
            # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
            p.val.Pathology <-
              max(lmodel$coefficients$mean[2, 4], 1e-323) #p.val.Pathology  #should not be zero
            coef.Pathology <- lmodel$coefficients$mean[2, 1]
            se.Pathology <- lmodel$coefficients$mean[2, 2]
            estimate.Pathology <- coef.Pathology #estimate.Pathology
            std.error.Pathology <- se.Pathology  #std.error.Pathology
            Position <- rownames(betas.j)
            
            return(c(
              p.val.Pathology,
              estimate.Pathology,
              std.error.Pathology,
              Position
            ))
          }
          
        }
        
      }
    
    ###Run EWAS using foreach() and %dopar% to tell R to run the analysis in parallel
    min.meth <- min(betas[betas > 0], na.rm = TRUE)
    max.meth <- max(betas[betas < 1], na.rm = TRUE)
    
    res <-
      foreach(
        j = 1:nrow(betas),
        .combine = rbind,
        .packages = c("betareg", "lmtest")
      ) %dopar% {
        testCpG(
          betas[j, ],
          pheno,
          formula = formula,
          link = link,
          max.meth = max.meth,
          min.meth = min.meth
        )
      }
    
    stopCluster(cl)
    
    colnames(res) <-
      c("p.val.Pathology",
        "estimate.Pathology",
        "std.error.Pathology",
        "Position")
    
    rownames(res) <- rownames(betas)
    
    res <- as.data.frame(res)
    head(res)
    
    return(res)
  }