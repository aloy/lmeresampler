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
  mod.rand <- model$call$random
  
  if(is.null(new.data)){
    # Place ystars in data
    mod.data <- model$data
    mod.data[,as.character(mod.fixd[[2]])] <- unname(new.y)
  } else{
    mod.data <- new.data
  }
  
  # create new lme
  ctrl <- nlme::lmeControl(opt = 'optim', returnObject = TRUE)
  if(is.null(mod.rand)){
    out.lme <- nlme::lme(fixed = mod.fixd, data = mod.data, control = ctrl)
  } else{
    mod.rand <- as.formula(mod.rand)
    out.lme <- nlme::lme(fixed = mod.fixd, data = mod.data, random = mod.rand, control = ctrl)
  }
  
  out.lme
}


#' Create list of Z matrices, similar to Ztlist in lme4
#' @importFrom forcats fct_inorder
#' @keywords internal
#' @noRd
extract_zlist.lme <- function(model){
  level.num <- ncol(model$groups)
  re.form <- formula(model$modelStruct$reStr)
  Z <- purrr::map(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  
  grp <- purrr::map(model$groups, forcats::fct_inorder)
  
  # if(level.num == 1) {
  #   Z <- as.data.frame(Z[[1]])
  #   Zlist <- purrr::map(Z, function(col) split(col, grp[[1]]))
  #   
  # } else{
    
    Z <- purrr::map(Z, as.data.frame)
    Z  <- Z[rev(names(Z))] # agree w/ order of model$group and bstar
    
    Zlist <- purrr::map(1:length(Z), function(i) purrr::map(Z[[i]], function(col) split(col, model$group[,i])))
    names(Zlist) <- names(Z)
  # }
  Zlist
}


.Zbstar.combine.lme <- function(bstar, Zlist){
  zbstar_list <- purrr::map(1:length(Zlist), function(e) {
    z.e <- Zlist[[e]]
    b.e <- bstar[[e]]
    purrr::map(1:length(z.e), function(j) unlist(mapply("*", z.e[[j]], b.e[,j], SIMPLIFY = FALSE)))
  })
  Reduce("+", unlist(zbstar_list, recursive = FALSE))
}


extract_parameters.lme <- function(model) {
  sig.e <- sigma(model)
  vc <- getVarCov(model)
  
  c(
    beta = fixef(model), 
    vc = c(diag(vc), sig.e^2)
  )
}
