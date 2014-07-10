#' Bootstrapping Linear Mixed Effects Models
#' 
#' \tabular{ll}{
#' Package: \tab lmeresampler\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.0\cr
#' Date: \tab 7/5/2014\cr
#' License: \tab GPLv3\cr
#' }
#' 
#' This is a package to help with bootstrapping Linear Mixed Effects Models.
#' 
#' @name lmeresampler
#' @docType package
#' @author Adam Loy and Spenser Steele \email{steeles@lawrence.edu}

library(lme4)
library(nlme)
library(roxygen)

#' @title Bootstrap for LMEs
#' 
#' @description
#' \code{bootstrap} helps streamline the bootstrap process for the parametric,
#' residual, cases, CGR, and REB bootstraps.
#' 
#' @details
#' 
#' @export
#' @param model The model to use
#' @param fn The function the user is interested in
#' @param type The \code{type} of bootstrap requested.
bootstrap <- function (model, fn, type, B){
  result <- switch(type,
         par = parametric.lmerMod(model, fn, B),
         res = residual.lmerMod(model, fn, B),
         case = case(model, fn, B),
         cgr = cgr(model, fn, B),
         reb = reb(model, fn, B, reb_type = 0),
         reb1 = reb(model, fn, B, reb_type = 1),
         reb2 = reb(model, fn, B, reb_type = 2))
  # TODO: need to be able to save results
}

#' @title Parametric Bootstrap
#' 
#' @description
#' The Parametric Bootstrap is uses the parametrically estimated
#' distribution function of the data to generate bootstrap samples.
#' 
#' @details
#' 
#' @inheritParams model
#' @inheritParams fn
#' @inheritParams B
#' 
#' @return list
parametric.lmerMod <- function(model, fn, B){
  fn <- match.fun(fn)
	
  model.fixef <- fixef(model) # Extract fixed effects
  y.star <- simulate(model, nsim = B)
  # Below is one idea that will be compatible with the boot package (for CIs)
  t0 <- fn(model)
  
  # Refit the model and apply 'fn' to it using lapply
  t.star <- lapply(y.star, function(x) {
    fn(refit(x, model))
  })
  
  t.star <- do.call("cbind", t.star) # Can these be nested?
  rownames(t.star) <- names(t0)
 
  RES <- structure(list(t0 = t0, t = t(t.star), R = B, data = model@frame, 
                        seed = .Random.seed, statistic = fn, 
                        sim = "parametric", call = match.call()), 
                        class = "boot")
  
  return(RES)
  
  # TODO: once we have things working, think about parallelization.
  #       using an llply statement would make this easy with the .parallel 
  #       parameter, but it might be slower than using mclapply, which is 
  #       found in the parallel package.
}


residual.lmerMod <- function (model, fn, B){
  fn <- match.fun(fn)
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # This needs to be run for every level
  level.num <- getME(object = model, name = "n_rfacs")

  
bstar <- lapply(model.ranef,
  FUN = function(x) {
    J <- nrow(x)
    
    # Sample of b*
    bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
    bstar <- x[bstar.index,]
    return(bstar)
  })

if(level.num == 1){
	bstar <- sapply(bstar, FUN = function(x) as.list(x))	
} 

  
#   for(i in 1:level.num){
#     temp.bstar <- calc_bstar(i)
#     bstar.vector <- c(bstar.vector, as.vector(t(temp.bstar)))
#   }
  
  Z <- getME(object = model, name = "Ztlist")
  Zbstar <- combine.elements(bstar = bstar, zstar = Z)

  
  # Sample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  # Combine function?
}

combine.elements <- function(bstar, zstar){
  lapply(1:length(), function(i){
    t(zstar[i]) %*% bstar[i]
  })
}

case <- function (model, fn){

}

cgr <- function (model, fn){
  
}

reb <- function (model, fn, reb_type = 1){
  if(reb_type = 1){
    # Call reb1 here
  }
  # reb code here
  if(reb_type = 2){
    # Call reb2 here
  }
}

reb1 <- function (model, fn){
  
}

reb2 <- function (model, fn){
  
}