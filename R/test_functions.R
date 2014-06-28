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
boot <- function (model, fn, type){
  switch(type,
         par = parametric(model, fn),
         res = residual(model, fn),
         case = case(model, fn),
         cgr = cgr(model, fn),
         reb = reb(model, fn, reb_type = 0),
         reb1 = reb(model, fn, reb_type = 1),
         reb2 = reb(model, fn, reb_type = 2))
}

#' @inheritParams model
#' @inheritParams fn
parametric <- function (model, fn){
  B <- 100 # Should be large
  D <- # Number of groups at level-2
  N <- # Number of samples in D groups
  for(b in 1:B){
    # Generate indpt. level-2 errors for D groups from a normal
    sigma.u <- # Find the estimated sigma_u^2
    u.star <- rnorm(D,0,sigma.u)
    
    # Generate indpt. level-1 errors for n samples from D groups
    sigma.e <- # Find the estimated sigma_e^2
    e.star <- rnorm(N,0,sigma.e)
    
    # Iterate through the entire list and simulate y.b using the model
    y.b <- rep(0,N)
    
    # Fit the model
  }
}

residual <- function (model, fn){
  
}

case <- function (model, fn){
  
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