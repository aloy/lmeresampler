#' @rdname bootstrap
bootstrap.lme <- function (model, fn, type, B, extra_step, reb_type){
  switch(type,
         par = parametric_bootstrap.lme(model, fn, B),
         resid = resid_bootstrap.lme(model, fn, B),
         case = case_bootstrap.lme(model, fn, B, extra_step = FALSE),
         cgr = cgr_bootstrap.lme(model, fn, B),
         reb = reb_bootstrap.lme(model, fn, B, reb_type = 0))
}

#' @rdname parametric_bootstrap
parametric_bootstrap.lme <- function(model, fn, B){
  # Match function
  fn <- match.fun(fn)
  # Extract fixed effects
  model.fixef <- nlme::fixef(model)
  
  ystar <- nlmeU::simulateY(model, nsim = B)
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

#' @rdname resid_bootstrap
resid_bootstrap.lme <- function(model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.resids(model = model)) )
  return(ystar)
}

#' @rdname cgr_bootstrap
cgr.lme <- function(model, fn, B){
  fn <- match.fun(fn)
  B <- 10
  
  ystar <- as.data.frame( replicate(n = B, .resample.cgr(model = model)) )
  
  return(ystar)
}


#' @rdname case_bootstrap
case.lme <- function (model, fn, B, extra_step = FALSE){
  # TODO: put everything below into lapply to replicate
  .cases.resamp <- function (model, extra_step){
    # Draw sample of size J from level-2 units
    model.split <- split(x=model$data, f=model$data[3])
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
  rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(model = model, extra_step = extra_step))))
  
  
  .cases.completion <- function(model, data, B, fn){
    t0 <- fn(model)
    
    # Refit the model and apply 'fn' to it using lapply
    form <- model@call$formula
    reml <- isREML(model)
    t.res <- matrix(0, ncol = 2, nrow = B)
    for(i in 1:B){
      myin <- ystar[,i]
      model.update <- update(object = model, fixed = myin ~ .)
      t.res[i,] <- fn(model.update)
    }
    t.res <- t(t.res)
    tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
    
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

#####################
# Utility Functions #
#####################

# Bootstrap completion
.bootstrap.completion.lme <- function(model, ystar, B, fn){
  t0 <- fn(model)
  
  t.res <- matrix(0, ncol = 2, nrow = B)
  for(i in 1:B){
    myin <- ystar[,i]
    model.update <- update(object = model, fixed = myin ~ .)
    t.res[i,] <- fn(model.update)
  }
  t.res <- t(t.res)
  tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
  
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  return(RES)
}

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
  Z.lme <- extract.lmeDesign(model)$Z
  one.Z <- matrix(1, ncol = ncol(Z.lme)/2, nrow = nrow(Z.lme))
  two.Z <- matrix(2, ncol = ncol(Z.lme)/2, nrow = nrow(Z.lme))
  my.counter <- 1
  for(i in 1:ncol(Z.lme)){
    if(i%%2==0){
      two.Z[,my.counter] <- Z.lme[,i]
      my.counter <- my.counter+1
    }else{
      one.Z[,my.counter] <- Z.lme[,i]}
    
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
  
  Z <- .extractZ.lme(model)
  
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
    names(ustar) <- names(Z)
  } else {
    ustar <- lapply(ustar, FUN = function(x) as.data.frame(x))
    ustar <- do.call(c, ustar)
    names(ustar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Get e*
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  # Combine
  ystar <- as.numeric(Xbeta + Zbstar.sum + estar)
  return(ystar)
}
