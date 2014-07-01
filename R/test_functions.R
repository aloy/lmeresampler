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

#' @param model The model to use
#' @param fn The function
#' @param type The type of bootstrap requested
#' @param FUN what function is the user interested in
boot <- function (model, fn, FUN, type){
  switch(type,
         par = parametric(model, fn, FUN),
         res = residual(model, fn, FUN),
         case = case(model, fn, FUN),
         cgr = cgr(model, fn, FUN),
         reb = reb(model, fn, FUN, reb_type = 0),
         reb1 = reb(model, fn, FUN, reb_type = 1),
         reb2 = reb(model, fn, FUN, reb_type = 2))
}

#' @inheritParams model
#' @inheritParams fn
parametric <- function (model, fn){
  B <- 100 # Should be large
  J <- # Number of groups at level-2
  N <- # Number of samples in J groups
  for(b in 1:B){
    # Generate indpt. level-2 errors for J groups from a normal
    sigma.u <- # Find the estimated sigma_u^2
    u.star <- rnorm(J,0,sigma.u)
    
    # Generate indpt. level-1 errors for n samples from J groups
    sigma.e <- # Find the estimated sigma_e^2
    e.star <- rnorm(N,0,sigma.e)
    
    # Iterate through the entire list and simulate y.b using the model
    y.b <- rep(0,N)
    
    # Fit the model
  }
}

#' @inheritParams model
#' @inheritParams fn
#' @inheritParams FUN
parametric2 <- function(model, fn, FUN, B){
  model.fixef <- fixef(model) # Extract fixed effects
  fn.star <- rep(0, B)
  model.star <- c(1:B)
  # Can I just simulate multiple times and run refit once on the entire array
  # to just get a returned array
  y.star <- simulate(model, nsim = B)
  # TODO: change refit according to notes in email
  # model.star[b] <- refit(model, y.star)
  model.star <- lapply(y.star, refit, object = model)
  # TODO: How to piece everything back together and what to return?
  
  return(model.star) # maybe a good thing to return?
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