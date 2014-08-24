fn <- match.fun(fn)
# Extract random effects
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
  
  level.num <- getME(object = model, name = "n_rfacs")
  
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
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}

ystar <- as.data.frame( replicate(n = B, .resample.cgr(model = model)) )

return(.bootstrap.completion(model, ystar, B, fn))