scale_center_uhat <- function(x){
  # Calculate the scaling factor for the ranefs
  S <- (t(x) %*% x) / nrow(x)
  R <- bdiag(lme4::VarCorr(model))
  Ls <- chol(S, pivot = TRUE)
  Lr <- chol(R, pivot = TRUE)
  A <- t(Lr %*% solve(Ls))
  
  # Rescale ranefs so empirical variance is equal to estimated variance
  Uhat <- x %*% A
  
  # zero center
  Uhat <- data.frame(scale(Uhat, scale = FALSE))
  
  return(Uhat)
}

scale_center_e <- function(x, sigma) {
  estar <- sigma * x %*% ((t(x) %*% x) / length(x))^(-1/2)
  scale(estar, scale = FALSE)
}


extract_parameters.lmerMod <- function(model) {
  sig.e <- sigma(model)
  vc <- as.data.frame(lme4::VarCorr(model))
  
  c(
    beta = getME(model, "beta"), 
    vc = vc$vcov[is.na(vc$var2)]
  )
}
