reb.lmerMod <- function (model, fn, B, reb_type = 0){
  
  fn <- match.fun(fn)
  
  
  
  .reb.two <- function(x) {
    #POST
    
    OneB <- matrix(c(1), nrow = B, ncol = 1)
    
    # calculate the log(var) of each
    u.lvar <- log(var(u))
    e.lvar <- log(var(e))
    # combine into a Bx2 matrix
    Sstar <- matrix(c(u.lvar, e.lvar), nrow = B, ncol = 2)
    # average of Sstar
    Mstar <- matrix(c(mean(u.lvar) * OneB, mean(e.lvar) * OneB), nrow = B, ncol = 2)
    # sd of Sstar
    Dstar <- matrix(c(sd(u.lvar) * OneB, sd(e.lvar) * OneB), nrow = B, ncol = 2)
    
    # Calculate 2x2 cov matrix of Sstar
    Cstar <- cov(Sstar)
    
    # Lstar
    Lstar <- Mstar + ((Sstar - Mstar) * Cstar^(-1/2)) * Dstar
    
    # Step c on pg 457
  }
  
  
  ystar <- as.data.frame( replicate(n = B, .resample.reb(model = model, reb_type = reb_type)) )
  
  # how do i extract the two sigma_u and sigma_e values from these results once
  # fitted below?
  
  # This step needs to be done outside the bootstrap
  if(reb_type == 2){
     # Refit the model and apply 'fn' to it using lapply
    tstar <- lapply(ystar, function(x) {
      fn(refit(object = model, newresp = x))
    })
    
    tstar <- do.call("cbind", tstar) # Can these be nested?
    rownames(tstar) <- names(fn(model))
    
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
    
    for(i in 1:length(estar)){
      estar[i,1] <- estar[i,1]-(mean(estar[,1])/length(estar[,1]))
    }
    for(i in 1:length(Uhat)){
      Uhat[i,1] <- Uhat[i,1]--(mean(Uhat[,1])/length(Uhat[,1]))
    }
    
  }else{
    Uhat <- u
    estar <- e
  }
  
  Xbeta <- predict(model, re.form = NA)
  
  # resample uhats
  
  Uhat <- as.data.frame(as.matrix(Uhat))
  
  Uhat.list <- list(Uhat)
  
  level.num <- getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    Uhat.list <- lapply(Uhat.list, FUN = function(x) as.list(x))[[1]]
    names(Uhat.list) <- names(Z)
  } else {
    Uhat.list <- sapply(Uhat.list, FUN = function(x) as.list(x))
  }
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  # Resample Uhat
  J <- length(Uhat.list[[1]])
  Uhat.samp <- sample(x = Uhat.list[[1]], size = J, replace = TRUE)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = as.data.frame(Uhat.samp), zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  
  return(y.star)
}



library(lme4)
library(HLMdiag)


# Data set
data(sleepstudy)

# model
(model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
B <- 10
fn <- fixef
