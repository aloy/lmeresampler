reb.lmerMod <- function (model, fn, B, reb_type = 0){
  
  fn <- match.fun(fn)
  
  
  
  .reb.two <- function(ue.mat) {
    #POST
    
    OneB <- matrix(c(1), nrow = B, ncol = 1)
    
    # combine into a Bx2 matrix
    Sstar <- as.matrix(ue.mat) # NOT A REAL MATRIX
    # average of Sstar
    Mstar <- Sstar %*% OneB
    # sd of Sstar
    Dstar <- matrix(c(sd(u.lvar) * OneB, sd(e.lvar) * OneB), nrow = B, ncol = 2)
    
    # Calculate 2x2 cov matrix of Sstar
    Cstar <- cov(Sstar)
    
    # Lstar
    Lstar <- Mstar + ((Sstar - Mstar) * Cstar^(-1/2)) * Dstar
    
    # Step c on pg 457
    return(Lstar)
  }
  
  
  ystar <- as.data.frame( replicate(n = B, .resample.reb(model = model, reb_type = reb_type)) )
  
  # how do i extract the two sigma_u and sigma_e values from these results once
  # fitted below?
  
  # This step needs to be done outside the bootstrap
  if(reb_type == 2){
     # Refit the model and apply 'fn' to it using lapply
    tstar <- lapply(ystar[1,], function(x) {
      fn(refit(object = model, newresp = x))
    })
    
    tstar <- do.call("cbind", tstar) # Can these be nested?
    rownames(tstar) <- names(fn(model))
    #NOT WORKING
    u.vec <- as.numeric(t(ystar[2,]))
    e.vec <- as.numeric(t(ystar[3,]))
    ue.mat <- matrix(c(u.vec, e.vec), ncol = 2)
    rt.res <- .reb.two(ue.mat)
    
  } else{
    RES <- .bootstrap.completion(model, ystar, B, fn)
  }
  
  return(RES)
}

.resample.reb <- function(model, reb_type){
  # use HLMresid to extract marginal residuals
  model.mresid <- HLMresid(object = model, type = "EB", level = "marginal")
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Z")
  
  # level 2 resid
  u <- solve(t(Z) %*% Z) %*% t(Z) %*% model.mresid
  # level 1 resid
  e <- model.mresid - Z %*% u
  if(reb_type == 1){
    #PRE
    
    
    # Calculations
    
    S <- (t(u)%*%u)/length(u)
    R <- bdiag(VarCorr(model))
    Ls <- chol(S, pivot = TRUE)
    Lr <- chol(R, pivot = TRUE)
    A <- t(Lr%*%solve(Ls))
    
    Uhat <- u%*%A
    
    sigma <- sigma(model)
    estar <- sigma * e %*% ((t(e) %*% e) / length(e))^(-1/2)
    
    # center
    estar <- scale(estar, scale = FALSE) # faster than the for loop
    Uhat <- scale(Uhat, scale = FALSE) 
    
  } else{
    Uhat <- u
    estar <- e
  }
  
  Xbeta <- predict(model, re.form = NA)
  
  # resample uhats
  
  Uhat <- as.data.frame(Uhat)
  Uhat.list <- list(Uhat)
  
  level.num <- getME(object = model, name = "n_rfacs")
  
  # Extract Z design matrix
  Ztlist <- getME(object = model, name = "Ztlist")
  
  if(level.num == 1){
    Uhat.list <- lapply(Uhat.list, FUN = function(x) as.list(x))[[1]]
    names(Uhat.list) <- names(Ztlist)
  } else {
    Uhat.list <- sapply(Uhat.list, FUN = function(x) as.list(x))
  }
    
  # Resample Uhat
  J <- length(Uhat.list[[1]])
  ustar <- sample(x = Uhat.list[[1]], size = J, replace = TRUE)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = as.data.frame(ustar), zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  # this is going to be a crude workaround
  u.lvar <- log(var(ustar))
  e.lvar <- log(var(estar))
  
  test.return <- structure(list(ystar = y.star, u = u.lvar, e = e.lvar))
  
  return(test.return)
}



library(lme4)
library(HLMdiag)


# Data set
data(sleepstudy)

# model
(model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
B <- 10
fn <- fixef
