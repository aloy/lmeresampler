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