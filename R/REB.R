reb.lmerMod <- function (model, fn, B, reb_type = 0){
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
    
    # center and scale level 2
    # I thought we didnt need to do this? (with the CGR at least)
    # center and scale level 1
    sigma <- sigma(model)
    estar <- sigma*e*((t(e)%*%e)/length(e))^(-1/2)
  }
  
  # Extract Z design matrix
  Z <- getME(object = model, name = "Ztlist")
  
  # level 2 resid
  rhbar <- tcrossprod(Z)^(-1) * t(Z) * model.mresid
  # level 1 resid
  
  # average the level 2 marginal resids
  model.mresid.avg <- sum(model.mresid) / length(model.mresid)
  
  
  
  
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
}






# Data set
data(sleepstudy)

# model
(model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
B <- 10
