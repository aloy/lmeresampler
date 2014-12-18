# Extract the residual covariance matrix from an lme object
.extractR.lme <- function(lme.fit) {
  n <- length( getResponse(lme.fit) )
  if (length(lme.fit$group) > 1) {
    stop("not implemented for multiple levels of nesting")
  } 
  else{
    ugroups <- unique(lme.fit$groups[[1]])
    if (!is.null(lme.fit$modelStruct$corStruct)) {
      V <- Matrix( corMatrix(lme.fit$modelStruct$corStruct) )
    }
    else V <- Diagonal(n)
  }
  if (!is.null(lme.fit$modelStruct$varStruct)) 
    sds <- 1/varWeights(lme.fit$modelStruct$varStruct)
  else sds <- rep(1, n)
  sds <- lme.fit$sigma * sds
  cond.var <- t(V * sds) * sds
  
  return(cond.var / lme.fit$sigma^2)
}

# Extract the ranef covariance matrix from an lme object
.extractD.lme <- function(lme.fit) {
  mod.mats <- .extract.lmeDesign(lme.fit)
  D <- Matrix( mod.mats$Vr )
  return(D)
}

# Extract the Z matrix from a model
.extractZ.lme <- function(model){
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
  Z <- structure(list(one = one.Z, two = two.Z))
  return(Z)
}

# Refit the model
updated.model<- function(model, new.y){
  # Extract formulas and data
  mod.fixd <- as.formula(model$call$fixed)
  mod.rand <- as.formula(model$call$random)
  mod.data <- model$data
  # Place ystars in data
  mod.data[,as.character(mod.fixd[[2]])] <- unname(new.y)
  # create new lme
  out.lme <- lme(fixed = mod.fixd, data = mod.data, random = mod.rand)
  return(out.lme)
}