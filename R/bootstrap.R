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
#' @author Adam Loy \email{loya@lawrence.edu}
#' @author Spenser Steele \email{steeles@lawrence.edu}

library(lme4)
library(nlme)
#library(roxygen)

#' @title Bootstrap for LMEs
#'
#' @description
#' \code{bootstrap} helps streamline the bootstrap process for the parametric,
#' residual, cases, CGR, and REB bootstraps.
#'
#' @details To choose a bootstrap use the \code{type} parameter.
#' For parametric use \code{"par"}, residual use \code{"res"}, cases use \code{"case"},
#' CGR use \code{"cgr"}, and REB use \code{"reb"}. The REB bootstrap has two types
#' which defaults to 1 but can be chosen using \code{"reb", "reb1", "reb2"}.
#'
#' @export
#' @param model The original model to use
#' @param fn The function the user is interested in testing
#' @param type The \code{type} of bootstrap requested, see details for types
#' @param B The number of bootstrap simulations
bootstrap <- function (model, fn, type, B){
  switch(type,
         par = parametric.lmerMod(model, fn, B),
         res = residual.lmerMod(model, fn, B),
         case = case.lmerMod(model, fn, B, extra_step = FALSE),
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
#' This function extracts the fixed effects, simulates from the model, refits the model
#' and then returns the results in a list.
#'
#' @inheritParams model
#' @inheritParams fn
#' @inheritParams B
#'
#' @return list
#'
#' @references
#'   @cite Chambers:2013ba
#'   @cite vanderLeeden:208kv
parametric.lmerMod <- function(model, fn, B){
  fn <- match.fun(fn)

  model.fixef <- fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)

  return(.bootstrap.completion(model, ystar, B, fn))

  # TODO: once we have things working, think about parallelization.
  #       using an llply statement would make this easy with the .parallel
  #       parameter, but it might be slower than using mclapply, which is
  #       found in the parallel package.
}


#' @title Residual Bootstrap
#'
#' @description
#' The Residual Bootstrap uses residuals to generate bootstrap samples.
#'
#' @details
#' This function extracts the Xbetas, random effects, residuals, and Z
#' design matrix in order to resample the residuals and complete the
#' bootstrap process.
#'
#' @inheritParams model
#' @inheritParams fn
#' @inheritParams B
#'
#' @return list
#'
#' @references
#'   @cite vanderLeeden:208kv
residual.lmerMod <- function (model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.resids(model = model)) )

  return(.bootstrap.completion(model, ystar, B, fn))
}

#' @title Cases Bootstrap
#'
#' @description
#' The Cases Bootstrap samples entire cases to generate the bootstrap.
#'
#' @details
#'
#' @param extra_step add the extra step
#' @inheritParams model
#' @inheritParams fn
#' @inheritParams B
#'
#' @return list
#'
#' @references
#'   @cite vanderLeeden:208kv
case.lmerMod <- function (model, fn, B, extra_step = FALSE){
  # TODO: put everything below into lapply to replicate
  .cases.resamp <- function (model, extra_step){
  # Draw sample of size J from level-2 units
  model.split <- split(x=model@frame, f=model@flist)
  model.split.samp <- sample(x=model.split, size = length(model.split), replace = TRUE)
  # For each sample, draw a sample of the cases from the level-2 unit
  
  if(extra_step == TRUE){
    model.resamp <- lapply(model.split.samp,
                           FUN = function(x) {
                             J <- nrow(x)
                             
                             # Sample of level-2 rows
                             model.sub.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                             resampled <- x[model.sub.index,]
                             return(resampled)
                           })
    model.comb <- do.call('rbind', model.resamp)
  } else{ # else statement needs to be located here
    model.comb <- do.call('rbind', model.split.samp)
  }
  }
  
  # DEPRECATED rep.data <- as.data.frame( replicate(n = B, .cases.resamp(model = model, extra_step = extra_step)) )
  rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(model = fm1, extra_step = extra_step))))
  
  
  .cases.completion <- function(model, data, B, fn){
    t0 <- fn(model)
    
    # Refit the model and apply 'fn' to it using lapply
    form <- model@call$formula
    reml <- isREML(model)
    tstar <- lapply(data, function(x) {
      fn(lmer(formula = form, data = x, REML = reml)) 
    })
    
    tstar <- do.call("cbind", tstar) # Can these be nested?
    rownames(tstar) <- names(t0)
    
    RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                          seed = .Random.seed, statistic = fn,
                          sim = "parametric", call = match.call()),
                     class = "boot")
    
    return(RES)
  }
  
  # Plugin to .cases.completion due to small changes
  return(.cases.completion(model, rep.data, B, fn))
}

#' @title CGR Bootstrap
#'
#' @description
#' 
#'
#' @details
#'
#' @inheritParams model
#' @inheritParams fn
#' @inheritParams B
#'
#' @return list
#'
#' @references
#'   @cite Chambers:2013ba
cgr.lmerMod <- function (model, fn, B){
  fn <- match.fun(fn)
  # Extract random effects
  model.ranef <- ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # scale the level-1 residuals and random effects
  # use two lapplys to do this over model.resid
  
  # center the scaled residuals at zero
  # also use an lapply here to do this, unless there is a faster way
  
  
  # Sample Random Effects
  ustar <- sample(x = model.ranef, size = length(model.ranef), replace = TRUE)
  # Resample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  
  # use model
  # fit
  # repeat B times
}
#####################
# Utility Functions #
#####################

#' @title Zbstar combine
#'
#' @description
#' Combine \code{bstar} and \code{zstar} to create {Zbstar}.
#'
#' @details
#' This function combines \code{bstar} and \code{zstar} to create {Zbstar} using an lapply statement.
#'
#' @param bstar A list of matrices bstar
#' @param zstar A list of matrices zstar
#'
#' @return matrix
.Zbstar.combine <- function(bstar, zstar){
  lapply(1:length(bstar), function(i){
    t(zstar[[i]]) %*% bstar[[i]]
  })
}

#' @title Bootstrap Completion
#'
#' @description
#' Finishes the bootstrap process and makes the output readable.
#'
#' @details
#' This function is given \code{model, ystar, B, fn} and uses them to complete
#' the bootstrap process. They are then structured into a list for output and returned.
#'
#' @param ystar The ystar being passed in
#' @inheritParams model
#' @inheritParams B
#' @inheritParams fn
#'
#' @return list
.bootstrap.completion <- function(model, ystar, B, fn){
  t0 <- fn(model)

  # Refit the model and apply 'fn' to it using lapply
  tstar <- lapply(ystar, function(x) {
    fn(refit(object = model, newresp = x))
  })

  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)

  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")

  return(RES)
}


#' Resampling residuals from mixed models

.resample.resids <- function(model){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  bstar <- lapply(model.ranef,
                  FUN = function(x) {
                    J <- nrow(x)
                    
                    # Sample of b*
                    bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                    bstar <- x[bstar.index,]
                    return(bstar)
                  })
  
  level.num <- getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    bstar <- lapply(bstar, FUN = function(x) as.list(x))[[1]]
    names(bstar) <- names(Z)
  } else {
    bstar <- sapply(bstar, FUN = function(x) as.list(x))
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  
  # Resample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
  
}