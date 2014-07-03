#' Helps with linear mixed effects resampling.
#' 
#' \tabular{ll}{
#' Package: \tab lmeresampler\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.0\cr
#' Date: \tab 6/27/2014\cr
#' License: \tab GPLv3\cr
#' }
#' 
#' Description here.
#' 
#' @name lmeresampler
#' @docType package
#' @author Adam Loy and Spenser Steele \email{steeles@lawrence.edu}

library(lme4)
library(nlme)
library(roxygen)

#' Bootstrap for LMEs
#' 
#' This is a short description of what the function does.
#' 
#' @export
#' @param model The model to use
#' @param fn The function
#' @param type The type of bootstrap requested
#' @param FUN what function is the user interested in
bootstrap <- function (model, fn, type, B){
  switch(type,
         par = parametric(model, fn, B),
         res = residual(model, fn, B),
         case = case(model, fn, B),
         cgr = cgr(model, fn, B),
         reb = reb(model, fn, B, reb_type = 0),
         reb1 = reb(model, fn, B, reb_type = 1),
         reb2 = reb(model, fn, B, reb_type = 2))
  # TODO: determine which is better a switch, or type <- match.arg(type)
  # TODO: need to be able to save results
}

#' @inheritParams model
#' @inheritParams fn

##############
# DEPRECATED #
##############
# parametric <- function (model, fn){
#   B <- 100 # Should be large
#   J <- # Number of groups at level-2
#   N <- # Number of samples in J groups
#   for(b in 1:B){
#     # Generate indpt. level-2 errors for J groups from a normal
#     sigma.u <- # Find the estimated sigma_u^2
#     u.star <- rnorm(J,0,sigma.u)
#     
#     # Generate indpt. level-1 errors for n samples from J groups
#     sigma.e <- # Find the estimated sigma_e^2
#     e.star <- rnorm(N,0,sigma.e)
#     
#     # Iterate through the entire list and simulate y.b using the model
#     y.b <- rep(0,N)
#     
#     # Fit the model
#   }
# }

#' @inheritParams model
#' @inheritParams fn
#' @inheritParams FUN
parametric <- function(model, fn, B){
  # QUESTION: are fn and FUN both for functions? If so, let's pick one.
  fn <- match.fun(fn)
	
  model.fixef <- fixef(model) # Extract fixed effects
  fn.star <- rep(0, B)
  # model.star <- c(1:B) # we don't need to initialize this is we use apply statements
  # Can I just simulate multiple times and run refit once on the entire array
  # to just get a returned array
  y.star <- simulate(model, nsim = B)
  model.star <- lapply(y.star, refit, object = model)
  # TODO: evaluate FUN for each refitted model to extract desired component.
  function(x) {
    fn(refit(x, model))
  }
  # Below is one idea that will be compatible with the boot package (for CIs)
  t0 <- fn(model)
  
  # Consider: lapply(___, FUN = function(x){fn(refit( ))})
  # llapply may be faster than lapply for nonparallelized code
  t.star <- lapply(model.star, fn)
  
  # Best option
  t.star <- lapply(y.star, function(x) {
    fn(refit(x, model))
  })
  
  # New lapply statement
  t.star <- lapply(lapply(y.star, refit, object = model), fn)
  # Or this
  t.star <- lapply(y.star, function(x) {
    fn(lapply(x, refit, object = model))
  })
  # Maybe this
  t.star <- lapply(y.star, function(x) {
    lapply(lapply(x, refit, object = model), fn)
  })
  
  t.star <- do.call("cbind", t.star) # Can these be nested?
  rownames(t.star) <- names(t0)
  
  # QUESTION: Is it faster to use only one lapply statement? It should be... try it.
  
  RES <- structure(list(t0 = t0, t = t(t.star), R = B, data = model@frame, 
                        seed = .Random.seed, statistic = fn, 
                        sim = "parametric", call = match.call()), 
                        class = "boot")
  
  return(RES) # maybe a good thing to return?
  
  # TODO: once we have things working, think about parallelization.
  #       using an llply statement would make this easy with the .parallel 
  #       parameter, but it might be slower than using mclapply, which is 
  #       found in the parallel package.
}


residual <- function (model, fn){
  B <- 100 # Should be large
  J <- # Number of groups at level-2
  N <- # Number of samples in J groups
  for(b in 1:B){
    # Draw sample from random effects of size J from estimated random effects
    # of level-2 residuals
    
    # Draw J sample residuals of size N with replacement from the level-1 residuals
    
    # Generate bootstrap samples
    
    # Compute estimates for all parameters of the two-level model
  }
}

case <- function (model, fn){
  B <- 100 # Should be large
  J <- # Number of groups at level-2
  N <- # Number of samples in J groups
  for(b in 1:B){
    # Draw sample of size J from level-2 units
  }
}

cgr <- function (model, fn){
  
}

reb <- function (model, fn, reb_type){
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