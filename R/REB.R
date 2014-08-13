reb.lmerMod <- function (model, fn, B, reb_type = 0){
  
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.reb(model = model, reb_type = reb_type)) )
  
  
  # This step needs to be done outside the bootstrap
  if(reb_type == 2){
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
  
  return(.bootstrap.completion(model, ystar, B, fn))
}

.resample.reb <- function(model, reb_type){
  # use HLMresid to extract marginal residuals
  model.mresid <- HLMresid(object = model, type = "EB", level = "marginal")
  # extract random effects
  model.ranef <- ranef(model)
  u <- as.matrix(model.ranef[[1]])
  # extract level 1 resids
  model.resid <- resid(model)
  e <- model.resid
  if(reb_type == 1){
    #PRE
    
    
    # Calculations
    
    ## HERE
    S <- (t(u)%*%u)/length(u)
    R <- bdiag(VarCorr(model))
    Ls <- chol(S, pivot = TRUE)
    Lr <- chol(R, pivot = TRUE)
    A <- t(Lr%*%solve(Ls))
    
    Uhat <- u%*%A
    ## To here might not be necessary b/c only working with level 2?
    
    sigma <- sigma(model)
    estar <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  }else{
    Uhat <- 
    estar <- 
  }
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  # level 2 resid
  # ISS: Running into issue here
  rhbar <- (t(Z) %*% Z)^(-1) * t(Z) * model.mresid
  # level 1 resid
  
  # average the level 2 marginal resids
  model.mresid.avg <- sum(model.mresid) / length(model.mresid)
  
  Xbeta <- predict(model, re.form = NA)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = Uhat, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  
  return(y.star)
}




# Data set
data(sleepstudy)

# model
(model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
B <- 10
