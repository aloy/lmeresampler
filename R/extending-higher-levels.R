.resample.cgr <- function(model){
  model.ranef <- ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  

  # Higher levels
  Uhat.list <- lapply(seq_along(model.ranef),
                 FUN = function(i) {
                   u <- scale(model.ranef[[i]], scale = FALSE)
                   S <- (t(u) %*% u) / length(u)
                   
                   re.name <- names(model.ranef)[i]
                   R <- bdiag(VarCorr(model)[[names(model.ranef)[i]]])
                  
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
  sigma <- sigma(model)
  ehat <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  Xbeta <- predict(model, re.form = NA)
  
#   Uhat <- as.data.frame(as.matrix(Uhat))
#   Uhat.list <- list(Uhat)
  
  level.num <- getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    Uhat.list <- lapply(Uhat.list, FUN = function(x) as.list(x))[[1]]
    names(Uhat.list) <- names(Z)
  } else {
    Uhat.list <- sapply(Uhat.list, FUN = function(x) as.list(x))
  }
  
  # Resample Uhat
  ustar <- lapply(Uhat.list,
                  FUN = function(df) {
                    sample(x = df, size = length(df), replace = TRUE)
                    })
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # sample
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}