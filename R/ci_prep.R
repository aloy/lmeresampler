# bootstrap CIs:
if(class(model) == "lmerMod"){
  
  ## normal t
  
  ### these are not in the right order, will need to be returned separately
  
  con <- data.frame(lme4::confint.merMod(model))
  colnames(con)[1] <- "norm.t.lower"
  norm.t.lower <- con[1]
  colnames(con)[2] <- "norm.t.upper"
  norm.t.upper <- con[2]
  
  ## boot t
  
  ### This will need to be changed for nlme because accessing coefficients is different
  
  # this gets all of the means (estimates) for boot_t calculation
  model.fits <- lme4::getME(model, "beta")
  model.fits <- append(model.fits, lme4::getME(model, "sigma"))
  model.fits <- append(model.fits, lme4::getME(model, "theta"))
  
  # sd of of estimates for boot_t calculation
  out <- summary(model)
  model.sds <- out$coefficients[,2] # fixef
  
  ## sd for variance components, thank you Stack Overflow!!
  dd.ML <- lme4:::devfun2(model,useSc=TRUE,signames=FALSE)
  
  vv <- as.data.frame(VarCorr(model)) ## need ML estimates!
  pars <- vv[,"sdcor"]
  ## will need to be careful about order if using this for
  ## a random-slopes model ...
  
  library("numDeriv")
  hh1 <- hessian(dd.ML,pars)
  vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
  ranef.sds <- sqrt(diag(vv2))  ## get standard errors
  
  model.sds <- append(model.sds, ranef.sds)
  
  # table of estimates and sds for boot_t calculation
  t.stats <- cbind(model.fits, model.sds)
  row.names(t.stats) <- colnames(replicates)
  t.stats <- cbind(t.stats, rep.mean)
  
  library(dplyr)
  t.stats <- as.data.frame(t.stats)
  t.stats <- t.stats %>% # thank you for the formula, Andy!!!
    mutate(boot_t  = (rep.mean - model.fits)/(model.sds/sqrt(B))) %>%
    mutate(boot.t.lower = ((model.fits - quantile(boot_t, 0.975)) * model.sds/sqrt(B))) %>%
    mutate(boot.t.upper = ((model.fits - quantile(boot_t, 0.025)) * model.sds/sqrt(B)))
  
  boot.t <- t.stats %>%
    select(boot.t.lower, boot.t.upper)
  
  ## percentile t
  
  perc.t.lower <- apply(replicates, 2, function(x) {
    round(quantile(x, 0.025), 8)
  })
  
  perc.t.upper <- apply(replicates, 2, function(x) {
    round(quantile(x, 0.975), 8)
  })
} else if(class(model) == "lme"){
  ## normal t
  
  ### intervals returns the same est. as doing summary() so let's use it!
  
  con <- nlme::intervals(vcmodB, level = 0.95, which = "all")
  norm.t.lower <- con$fixed[, 1]  # lower interval for all fixedef
  norm.t.lower <- append(norm.t.lower, con$reStruct[[1]][1]) # lower interval for ranef
  norm.t.lower <- append(norm.t.lower, con$sigma[1]) # lower interval for sigma
  norm.t.lower <- unlist(norm.t.lower)
  names(t.stats) <- colnames(replicates) # something needs to be done about the names
  
  norm.t.upper <- con$fixed[, 3]  # upper interval for all fixedef
  norm.t.upper <- append(norm.t.upper, con$reStruct[[1]][3]) # upper interval for ranef
  norm.t.upper <- append(norm.t.upper, con$sigma[3]) # upper interval for sigma
  norm.t.upper <- unlist(norm.t.upper)
  names(norm.t.upper) <- colnames(replicates) # something needs to be done about the names
  
  ## boot t
  
  # this gets all of the means (estimates) for boot_t calculation
  model.fits <- con$fixed[, 2]  # estimates for all fixef
  model.fits <- append(model.fits, con$reStruct[[1]][2]) # estimate for ranef
  model.fits <- append(model.fits, con$sigma[2]) # estimate for sigma
  model.fits <- unlist(model.fits)
  names(model.fits) <- colnames(replicates) # something needs to be done about the names
  
  # sd of of estimates for boot_t calculation
  out <- summary(model)
  model.sds <- out$coefficients[,2] # fixef
  
  ## sd for variance components, thank you Stack Overflow!!
  dd.ML <- lme4:::devfun2(model,useSc=TRUE,signames=FALSE)
  
  vv <- as.data.frame(VarCorr(model)) ## need ML estimates!
  pars <- vv[,"sdcor"]
  ## will need to be careful about order if using this for
  ## a random-slopes model ...
  
  library("numDeriv")
  hh1 <- hessian(dd.ML,pars)
  vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
  ranef.sds <- sqrt(diag(vv2))  ## get standard errors
  
  model.sds <- append(model.sds, ranef.sds)
  
  # table of estimates and sds for boot_t calculation
  t.stats <- cbind(model.fits, model.sds)
  row.names(t.stats) <- colnames(replicates)
  t.stats <- cbind(t.stats, rep.mean)
  
  library(dplyr)
  t.stats <- as.data.frame(t.stats)
  t.stats <- t.stats %>% # thank you for the formula, Andy!!!
    mutate(boot_t  = (rep.mean - model.fits)/(model.sds/sqrt(B))) %>%
    mutate(boot.t.lower = ((model.fits - quantile(boot_t, 0.975)) * model.sds/sqrt(B))) %>%
    mutate(boot.t.upper = ((model.fits - quantile(boot_t, 0.025)) * model.sds/sqrt(B)))
  
  boot.t <- t.stats %>%
    select(boot.t.lower, boot.t.upper)
  
  ## percentile t
  
  perc.t.lower <- apply(replicates, 2, function(x) {
    round(quantile(x, 0.025), 8)
  })
  
  perc.t.upper <- apply(replicates, 2, function(x) {
    round(quantile(x, 0.975), 8)
  })
}

