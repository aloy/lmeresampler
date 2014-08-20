.resample.cgr <- function(model){
  model.ranef <- ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Level 2
  u <- model.ranef
  
  # Calculations
  S <- (t(u)%*%u)/length(u)
  R <- bdiag(VarCorr(model))
  Ls <- chol(S, pivot = TRUE)
  Lr <- chol(R, pivot = TRUE)
  A <- t(Lr%*%solve(Ls))
  
  Uhat <- as.matrix(u%*%A)
  Uhat <- as.data.frame(Uhat)
  
  
  # Level 1
  e <- model.resid
  sigma <- sigma(model)
  ehat <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  Xbeta <- predict(model, re.form = NA)
  
  Uhat <- as.data.frame(as.matrix(Uhat))
  Uhat.list <- list(Uhat)
  
  level.num <- getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    Uhat.list <- lapply(Uhat.list, FUN = function(x) as.list(x))[[1]]
    names(Uhat.list) <- names(Z)
  } else {
    Uhat.list <- sapply(Uhat.list, FUN = function(x) as.list(x))
  }
  
  # Resample Uhat
  ustar <- sample(x = Uhat.list[[1]], size = length(Uhat.list[[1]]), replace = TRUE)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = as.data.frame(ustar), zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # sample
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}