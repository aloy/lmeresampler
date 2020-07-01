# Extract the residual covariance matrix from an lme object
.extractR.lme <- function(lme.fit) {
  n <- length( nlme::getResponse(lme.fit) )
  if (length(lme.fit$group) > 1) {
    stop("not implemented for multiple levels of nesting")
  } 
  else{
    ugroups <- unique(lme.fit$groups[[1]])
    if (!is.null(lme.fit$modelStruct$corStruct)) {
      V <- Matrix( nlme::corMatrix(lme.fit$modelStruct$corStruct) )
    }
    else V <- Diagonal(n)
  }
  if (!is.null(lme.fit$modelStruct$varStruct)) 
    sds <- 1/nlme::varWeights(lme.fit$modelStruct$varStruct)
  else sds <- rep(1, n)
  sds <- lme.fit$sigma * sds
  cond.var <- t(V * sds) * sds
  
  return(cond.var / lme.fit$sigma^2)
}

# Extract the ranef covariance matrix from an lme object
.extractD.lme <- function(lme.fit) {
  mod.mats <- RLRsim::extract.lmeDesign(lme.fit)
  D <- Matrix( mod.mats$Vr )
  return(D)
}

# Extract the Z matrix from a model
.extractZ.lme <- function(model){
  Z.lme <- RLRsim::extract.lmeDesign(model)$Z
  one.Z <- matrix(1, ncol = ncol(Z.lme)/2, nrow = nrow(Z.lme))
  two.Z <- matrix(2, ncol = ncol(Z.lme)/2, nrow = nrow(Z.lme))
  my.counter <- 1
  for(i in 1:ncol(Z.lme)){
    if(i%%2==0){
      two.Z[,my.counter] <- Z.lme[,i]
      my.counter <- my.counter+1
    }else{
      one.Z[,my.counter] <- Z.lme[,i]}
    
  }
  one.Z <- t(one.Z)
  two.Z <- t(two.Z)
  Z <- structure(list(one = one.Z, two = two.Z))
  return(Z)
}

# Refit the model
updated.model<- function(model, new.y = NULL, new.data = NULL){
  # Extract formulas and data
  mod.fixd <- as.formula(model$call$fixed)
  mod.rand <- as.formula(model$call$random)
  
  if(is.null(new.data)){
    # Place ystars in data
    mod.data <- model$data
    mod.data[,as.character(mod.fixd[[2]])] <- unname(new.y)
  } else{
    mod.data <- new.data
  }
  
  # create new lme
  ctrl <- nlme::lmeControl(opt = 'optim')
  out.lme <- nlme::lme(fixed = mod.fixd, data = mod.data, random = mod.rand, control = ctrl)
  return(out.lme)
}

