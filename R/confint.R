#' @rdname confint
#' @export
# bootstrap CI method for object of class lmeresamp
confint.lmeresamp <- function(object, method, level) {
  
  if(missing(level)){
    level <- 0.95
  } else if(!level %in% (0:1)){
    stop("please specify a confidence level between 0 and 1")
  }
  
  if(class(object$model) == "lmerMod"){ ## normal t
    
    if(method == "norm") {
      
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      colnames(con)[1] <- "norm.t.lower"
      norm.t.lower <- con[1]
      colnames(con)[2] <- "norm.t.upper"
      norm.t.upper <- con[2]
      
      .norm.t.completion(norm.t.lower, norm.t.upper)
      
    } else if(method == "boot-t"){ ## boot t
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- lme4::getME(object$model, "beta")
      model.fits <- append(model.fits, lme4::getME(object$model, "sigma"))
      model.fits <- append(model.fits, lme4::getME(object$model, "theta"))
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$coefficients[, 2] # fixef
      
      ### sd for variance components, thank you Ben Bolker!!
      dd.ML <- lme4:::devfun2(object$model, useSc=TRUE, signames=FALSE)
      
      vv <- as.data.frame(VarCorr(object$model)) ## need ML estimates!
      pars <- vv[,"sdcor"]
      ## will need to be careful about order if using this for
      ## a random-slopes model ...
      
      library("numDeriv")
      hh1 <- hessian(dd.ML,pars)
      vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
      ranef.sds <- sqrt(diag(vv2))  ## get standard errors
      
      model.sds <- append(model.sds, ranef.sds)
      
      .boot.t.completion(object, level, model.fits, model.sds)
      
    } else if(method == "perc"){ ## percentile t
      
      .perc.t.completion(object, level)
      
    } else if(method == "all"){
      
      ## normal-t
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      colnames(con)[1] <- "norm.t.lower"
      norm.t.lower <- con[1]
      colnames(con)[2] <- "norm.t.upper"
      norm.t.upper <- con[2]
      
      .norm.t.completion(norm.t.lower, norm.t.upper)
      
      ## boot-t
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- lme4::getME(object$model, "beta")
      model.fits <- append(model.fits, lme4::getME(object$model, "sigma"))
      model.fits <- append(model.fits, lme4::getME(object$model, "theta"))
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$coefficients[, 2] # fixef
      
      ### sd for variance components, thank you Ben Bolker!!
      dd.ML <- lme4:::devfun2(object$model, useSc=TRUE, signames=FALSE)
      
      vv <- as.data.frame(VarCorr(object$model)) ## need ML estimates!
      pars <- vv[,"sdcor"]
      ## will need to be careful about order if using this for
      ## a random-slopes model ...
      
      hh1 <- hessian(dd.ML,pars)
      vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
      ranef.sds <- sqrt(diag(vv2))  ## get standard errors
      
      model.sds <- append(model.sds, ranef.sds)
      
      .boot.t.completion(object, level, model.fits, model.sds)
      
      ## percentile-t
      .perc.t.completion(object, level)
      
    } else{
      stop("'method' must be either 'norm', 'boot-t', 'perc', or 'all'")
    }
  } else if(class(model) == "lme"){
    
    if(method == "norm"){ ## normal t
      
      ### intervals returns the same estimates as doing summary() so let's use it!
      con <- nlme::intervals(object$model, level = level, which = "all")
      norm.t.lower <- con$fixed[, 1]  # lower interval for all fixedef
      norm.t.lower <- append(norm.t.lower, con$reStruct[[1]][1]) # lower interval for ranef
      norm.t.lower <- append(norm.t.lower, con$sigma[1]) # lower interval for sigma
      norm.t.lower <- unlist(norm.t.lower)
      names(norm.t.lower)[(length(con$fixed[, 1]) + 1) : (length(con$fixed[, 1]) + length(con$reStruct[[1]][1]))] <- row.names(con$reStruct[[1]][1])
      names(norm.t.lower)[(length(con$fixed[, 1]) + length(con$reStruct[[1]][1]) + 1)] <- "sigma"
      
      norm.t.upper <- con$fixed[, 3]  # upper interval for all fixedef
      norm.t.upper <- append(norm.t.upper, con$reStruct[[1]][3]) # upper interval for ranef
      norm.t.upper <- append(norm.t.upper, con$sigma[3]) # upper interval for sigma
      norm.t.upper <- unlist(norm.t.upper)
      names(norm.t.upper)[(length(con$fixed[, 3]) + 1) : (length(con$fixed[, 3]) + length(con$reStruct[[1]][3]))] <- row.names(con$reStruct[[1]][3])
      names(norm.t.upper)[(length(con$fixed[, 3]) + length(con$reStruct[[1]][3]) + 1)] <- "sigma"
      
      .norm.t.completion(norm.t.lower, norm.t.upper)
      
    } else if(method == "boot-t"){  ## boot t
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- con$fixed[, 2]  # estimates for all fixef
      model.fits <- append(model.fits, con$reStruct[[1]][2]) # estimate for ranef
      model.fits <- append(model.fits, con$sigma[2]) # estimate for sigma
      model.fits <- unlist(model.fits)
      names(model.fits)[(length(con$fixed[, 2]) + 1) : (length(con$fixed[, 2]) + length(con$reStruct[[1]][2]))] <- row.names(con$reStruct[[1]][2])
      names(model.fits)[(length(con$fixed[, 2]) + length(con$reStruct[[1]][2]) + 1)] <- "sigma"
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$tTable[, 2] # fixef 
      
      ## construct deviance function
      devfun <- do.call(mkLmerDevfun, object$model)
      
      ### sd for variance components, thank you Ben Bolker!!
      getTheta <- function(phi,sigma,nmax) {
        ## make corStruct: fake data sequence within a single block
        cc <- nlme::Initialize(nlme::corAR1(phi),data=data.frame(t=seq(nmax)))
        ## get inverse Cholesky factor
        mm <- matrix(nlme::corFactor(cc),nrow=nmax) ## 
        ## check backsolve() idiom: all.equal(solve(mm),backsolve(mm,diag(nmax),upper.tri=FALSE))
        mm2 <- backsolve(mm,diag(nmax),upper.tri=FALSE) ## invert ...
        return(sigma*mm2[lower.tri(mm2,diag=TRUE)])     ## take lower tri & scale by SD
      }
      
      devfun2.2 <- function(theta,nmax) { # no idea what nmax is
        new_theta <- getTheta(phi=theta[2],sigma=theta[3],nmax)
        devfun(c(theta[1],new_theta))
      }
      
      dd.ML <- devfun2.2(c(1,0.5,1),nmax=20)
      
      vv <- VarCorr(object$model) ## need ML estimates!
      pars <- as.numeric(vv[,"StdDev"]) # not sure if losing the label names is bad
      ## will need to be careful about order if using this for
      ## a random-slopes model ...
      
      hh1 <- hessian(dd.ML,pars)
      vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
      ranef.sds <- sqrt(diag(vv2))  ## get standard errors
      
      model.sds <- append(model.sds, ranef.sds)
      
      .boot.t.completion(object, level, model.fits, model.sds)
      
    } else if(method == "perc"){ ## percentile t
      
      .perc.t.completion(object, level)
      
    }
  } else if(method == "all"){
    
    ## normal-t
    ### intervals returns the same estimates as doing summary() so let's use it!
    con <- nlme::intervals(object$model, level = level, which = "all")
    norm.t.lower <- con$fixed[, 1]  # lower interval for all fixedef
    norm.t.lower <- append(norm.t.lower, con$reStruct[[1]][1]) # lower interval for ranef
    norm.t.lower <- append(norm.t.lower, con$sigma[1]) # lower interval for sigma
    norm.t.lower <- unlist(norm.t.lower)
    names(norm.t.lower)[(length(con$fixed[, 1]) + 1) : (length(con$fixed[, 1]) + length(con$reStruct[[1]][1]))] <- row.names(con$reStruct[[1]][1])
    names(norm.t.lower)[(length(con$fixed[, 1]) + length(con$reStruct[[1]][1]) + 1)] <- "sigma"
    
    norm.t.upper <- con$fixed[, 3]  # upper interval for all fixedef
    norm.t.upper <- append(norm.t.upper, con$reStruct[[1]][3]) # upper interval for ranef
    norm.t.upper <- append(norm.t.upper, con$sigma[3]) # upper interval for sigma
    norm.t.upper <- unlist(norm.t.upper)
    names(norm.t.upper)[(length(con$fixed[, 3]) + 1) : (length(con$fixed[, 3]) + length(con$reStruct[[1]][3]))] <- row.names(con$reStruct[[1]][3])
    names(norm.t.upper)[(length(con$fixed[, 3]) + length(con$reStruct[[1]][3]) + 1)] <- "sigma"
    
    .norm.t.completion(norm.t.lower, norm.t.upper)
    
    ## boot-t
    ### this gets all of the means (estimates) for boot_t calculation
    model.fits <- con$fixed[, 2]  # estimates for all fixef
    model.fits <- append(model.fits, con$reStruct[[1]][2]) # estimate for ranef
    model.fits <- append(model.fits, con$sigma[2]) # estimate for sigma
    model.fits <- unlist(model.fits)
    names(model.fits)[(length(con$fixed[, 2]) + 1) : (length(con$fixed[, 2]) + length(con$reStruct[[1]][2]))] <- row.names(con$reStruct[[1]][2])
    names(model.fits)[(length(con$fixed[, 2]) + length(con$reStruct[[1]][2]) + 1)] <- "sigma"
    
    ### sd of of estimates for boot_t calculation
    out <- summary(object$model)
    model.sds <- out$tTable[,2] # fixef 
    
    ### sd for variance components, thank you Ben Bolker!!
    getTheta <- function(phi,sigma,nmax) {
      ## make corStruct: fake data sequence within a single block
      cc <- nlme::Initialize(nlme::corAR1(phi),data=data.frame(t=seq(nmax)))
      ## get inverse Cholesky factor
      mm <- matrix(nlme::corFactor(cc),nrow=nmax) ## 
      ## check backsolve() idiom: all.equal(solve(mm),backsolve(mm,diag(nmax),upper.tri=FALSE))
      mm2 <- backsolve(mm,diag(nmax),upper.tri=FALSE) ## invert ...
      return(sigma*mm2[lower.tri(mm2,diag=TRUE)])     ## take lower tri & scale by SD
    }
    
    devfun2.2 <- function(theta,nmax) { # no idea what nmax is
      new_theta <- getTheta(phi=theta[2],sigma=theta[3],nmax)
      devfun(c(theta[1],new_theta))
    }
    
    dd.ML <- devfun2.2(object$model, useSc=TRUE, signames=FALSE)
    
    vv <- VarCorr(object$model) ## need ML estimates!
    pars <- as.numeric(vv[,"StdDev"]) # not sure if losing the label names is bad
    ## will need to be careful about order if using this for
    ## a random-slopes model ...
    
    hh1 <- hessian(dd.ML,pars)
    vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
    ranef.sds <- sqrt(diag(vv2))  ## get standard errors
    
    model.sds <- append(model.sds, ranef.sds)
    
    .boot.t.completion(object, level, model.fits, model.sds)
    
    ## percentile t
    .perc.t.completion(object, level)
    
  } else{
    stop("'method' must be either 'norm', 'boot-t', 'perc', or 'all'")
  }
}


#' @title Percentile-t interval completion
#'
#' @description
#' Execute the percentile-t interval process
#'
#' @details
#' This function uses \code{object} and \code{level} to calculate a percentile-t
#' interval for the fixed and random components
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.perc.t.completion <- function(object, level){
  
  perc.t.lower <- apply(object$replicates, 2, function(x) {
    round(quantile(x, (1 - level)/2), 8)
  })
  
  perc.t.upper <- apply(object$replicates, 2, function(x) {
    round(quantile(x, level + (1 - level)/2), 8)
  })
  
  perc.t <- data.frame(cbind(perc.t.lower, perc.t.upper))
  cat(paste("95% percentile-t interval: \n"))
  print(perc.t)
  cat(paste("\n"))
}

#' @title Bootstrap-t interval completion
#'
#' @description
#' Finish the bootstrap-t interval process
#'
#' @details
#' This function uses \code{object} and \code{level} to calculate a bootstrap-t
#' interval for the fixed and random components
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.boot.t.completion <- function(object, level, model.fits, model.sds){
  
  ### table of estimates and sds for boot_t calculation
  t.stats <- cbind(model.fits, model.sds)
  
  row.names(t.stats) <- colnames(object$replicates)
  t.stats <- cbind(t.stats, object$stats$rep.mean)
  
  t.stats <- as.data.frame(t.stats)
  t.stats <- t.stats %>% # thank you for the formula, Andy!!!
    mutate(boot_t  = (object$stats$rep.mean - model.fits)/(model.sds/sqrt(object$R))) %>%
    mutate(boot.t.lower = ((model.fits - quantile(boot_t, level + (1 - level)/2)) * model.sds/sqrt(object$R))) %>%
    mutate(boot.t.upper = ((model.fits - quantile(boot_t, (1 - level)/2)) * model.sds/sqrt(object$R)))
  
  boot.t <- t.stats %>%
    select(boot.t.lower, boot.t.upper)
  
  cat(paste("95% bootstrap-t interval: \n"))
  print(boot.t)
  cat(paste("\n"))
}

#' @title Normal-t interval completion
#'
#' @description
#' Complete the normal-t interval process
#'
#' @details
#' This function uses \code{norm.t.lower} and \code{norm.t.upper} to make a normal-t
#' interval for the fixed and random components
#'
#' @param norm.t.lower The lower bound of the interval
#' @param norm.t.upper The upper bound of the interval
#'
#' @keywords internal
#' @noRd
.norm.t.completion <- function(norm.t.lower, norm.t.upper){
  
  norm.t <-  data.frame(cbind(norm.t.lower, norm.t.upper))
  
  cat(paste("95% normal-t interval: \n"))
  print(norm.t)
  cat(paste("\n"))
}
