fn <- match.fun(fn)
B <- 10
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
                        vc.temp <- getVarCov(model)
                        R <- matrix(0, ncol = ncol(vc.temp), nrow = nrow(vc.temp))
#                         for(i in 1:ncol(vc.temp)){
#                           for(j in 1:nrow(vc.temp)){
#                             R[j,i] <- vc.temp[j,i]
#                           }
#                         }
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
  sigma <- sigma(model)
  ehat <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  
  # Extract and construct Z design matrix
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
  Z.str <- structure(list(one = one.Z, two = two.Z))
  
  Xbeta <- predict(model, re.form = NA)
  
  #level.num <- getME(object = model, name = "n_rfacs")
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
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}

ystar <- as.data.frame( replicate(n = B, .resample.cgr(model = model)) )