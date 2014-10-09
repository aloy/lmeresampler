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

#library(lme4)
library(nlme)
library(HLMdiag)
library(RLRsim)
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
bootstrap.nlme <- function (model, fn, type, B){
  switch(type,
         par = parametric.nlme(model, fn, B),
         res = residual.nlme(model, fn, B),
         case = case.nlme(model, fn, B, extra_step = FALSE),
         cgr = cgr.nlme(model, fn, B),
         reb = reb.nlme(model, fn, B, reb_type = 0))
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
parametric.nlme <- function(model, fn, B){
  # Match function
  fn <- match.fun(fn)
  # Extract fixed effects
  model.fixef <- fixed.effects(model)
  
  ystar <- simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  
  t0 <- fn(model)
  
  t.res <- list()
  length(t.res) <- B
  #    t.res <- matrix(0, ncol = 2, nrow = B)
  #   for(i in 1:B){
  #     myin <- ystar[,i]
  #     model.update <- update(object = model, fixed = myin ~ .)
  #     t.res[i,] <- fn(model.update)
  #   }
  
  for(i in 1:B){
    #     myin <- ystar[,i]
    #     model.update <- nlme:::update.lme(object = model, fixed = myin ~ .)
    #     t.res[i,] <- fn(model.update)
    try.fit <- try( updated.model(model = model, new.y = ystar[,i]) )
    ### NOTE: Need to check if there was an error in try, and then use fn()
    t.res[[i]] <- updated.model(model = model, new.y = ystar[,i])
    
  }
  t.res <- do.call('rbind', t.res)
  tstar <- data.frame(t(t.res))
  
  #   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
  #   
  #   tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  return(RES)
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
residual.nlme <- function(model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.resids(model = model)) )
  return(ystar)
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
cgr.nlme <- function(model, fn, B){
  fn <- match.fun(fn)
  B <- 10
  
  ystar <- as.data.frame( replicate(n = B, .resample.cgr(model = model)) )
  
  return(ystar)
}

#####################
# Utility Functions #
#####################

# Extract the residual covariance matrix from an lme object
.extractR.lme <- function(lme.fit) {
  n <- length( getResponse(lme.fit) )
  if (length(lme.fit$group) > 1) {
    stop("not implemented for multiple levels of nesting")
  } 
  else{
    ugroups <- unique(lme.fit$groups[[1]])
    if (!is.null(lme.fit$modelStruct$corStruct)) {
      V <- Matrix( corMatrix(lme.fit$modelStruct$corStruct) )
    }
    else V <- Diagonal(n)
  }
  if (!is.null(lme.fit$modelStruct$varStruct)) 
    sds <- 1/varWeights(lme.fit$modelStruct$varStruct)
  else sds <- rep(1, n)
  sds <- lme.fit$sigma * sds
  cond.var <- t(V * sds) * sds
  
  return(cond.var / lme.fit$sigma^2)
}

# Extract the ranef covariance matrix from an lme object
.extractD.lme <- function(lme.fit) {
  mod.mats <- .extract.lmeDesign(lme.fit)
  D <- Matrix( mod.mats$Vr )
  return(D)
}

# Extract the Z matrix from a model
.extractZ.lme <- function(model){
  Z.nlme <- extract.lmeDesign(model)$Z
  one.Z <- matrix(1, ncol = ncol(Z.nlme)/2, nrow = nrow(Z.nlme))
  two.Z <- matrix(2, ncol = ncol(Z.nlme)/2, nrow = nrow(Z.nlme))
  my.counter <- 1
  for(i in 1:ncol(Z.nlme)){
    if(i%%2==0){
      two.Z[,my.counter] <- Z.nlme[,i]
      my.counter <- my.counter+1
    }else{
      one.Z[,my.counter] <- Z.nlme[,i]}
    
  }
  one.Z <- t(one.Z)
  two.Z <- t(two.Z)
  Z <- structure(list(one = one.Z, two = two.Z))
  return(Z)
}

# Refit the model
updated.model<- function(model, new.y){
  # Extract formulas and data
  mod.fixd <- as.formula(model$call$fixed)
  mod.rand <- as.formula(model$call$random)
  mod.data <- model$data
  # Place ystars in data
  mod.data[,as.character(mod.fixd[[2]])] <- unname(new.y)
  # create new lme
  out.lme <- lme(fixed = mod.fixd, data = mod.data, random = mod.rand)
  return(out.lme)
}

#' Resampling residuals from mixed models
#'
#' @inheritParams model
#' 
.resample.resids <- function(model){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- random.effects(model)
  
  # Extract residuals
  model.resid <- residuals(model)
  
  # Extract Z design matrix
  # Z <- getME(object = model, name = "Ztlist") # FIX THIS
  
  J <- nrow(model.ranef)
  
  # Create Z design matrix
  Z <- diag(J)
  
  # Sample of b*
  bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
  bstar <- model.ranef[bstar.index,]
  
  level.num <- model$dims$Q # I believe this is where it lists it?
  
  if(level.num == 1){
    bstar <- lapply(bstar, FUN = function(x) as.list(x))[[1]]
    #names(bstar) <- names(Z)
  } else {
    bstar <- lapply(bstar, FUN = function(x) as.data.frame(x))
    bstar <- do.call(c, bstar)
    #names(bstar) <- names(Z)
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

#' CGR resampling procedures
#' 
#'
#' @inheritParams model
.resample.cgr <- function(model){
  model.ranef <- random.effects(model)
  
  # Extract residuals
  model.resid <- residuals(model)
  
  
  # Higher levels
  Uhat.list <- lapply(seq_along(model.ranef),
                      FUN = function(i) {
                        u <- scale(model.ranef[[i]], scale = FALSE)
                        S <- (t(u) %*% u) / length(u)
                        
                        re.name <- names(model.ranef)[i]
                        vc.temp <- getVarCov(model)
                        R <- getVarCov(model)[1]
                        
                        Ls <- chol(S, pivot = TRUE)
                        Lr <- chol(R, pivot = TRUE)
                        A <- t(Lr %*% solve(Ls))
                        
                        Uhat <- as.matrix(u %*% A)
                        Uhat <- as.data.frame(Uhat)
                        
                        return(Uhat)
                      })  
  names(Uhat.list) <- names(model.ranef)
  
  # Level 1
  e <- as.numeric(scale(model.resid, scale = FALSE))
  sigma <- model$sigma
  ehat <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  
  # Extract and construct Z design matrix
  
  Z.str <- .extractZ.lme(model)
  
  Xbeta <- predict(model, re.form = NA)
  
  level.num <- 1
  
  # Resample Uhat
  ustar <- lapply(Uhat.list,
                  FUN = function(df) {
                    index <- sample(x = seq_len(nrow(df)), size = nrow(df), replace = TRUE)
                    return(df[index,])
                  })
  
  # Structure u*
  if(level.num == 1){
    if(is.data.frame(ustar[[1]])){
      ustar <- lapply(ustar, FUN = function(x) as.list(x))[[1]] 
    }
    names(ustar) <- names(Z.str)
  } else {
    ustar <- lapply(ustar, FUN = function(x) as.data.frame(x))
    ustar <- do.call(c, ustar)
    names(ustar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Z.str)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Get e*
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  # Combine
  ystar <- as.numeric(Xbeta + Zbstar.sum + estar)
  return(ystar)
}
