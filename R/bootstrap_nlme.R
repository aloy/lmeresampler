#' @rdname bootstrap
#' @export
bootstrap.lme <- function (model, fn, type, B, resample, reb_type){
  switch(type,
         parametric = parametric_bootstrap.lme(model, fn, B),
         residual = resid_bootstrap.lme(model, fn, B),
         case = case_bootstrap.lme(model, fn, B, resample),
         cgr = cgr_bootstrap.lme(model, fn, B),
         reb = reb_bootstrap.lme(model, fn, B, reb_type = 0))
}


#' @rdname parametric_bootstrap
#' @export
#' @importFrom nlmeU simulateY
parametric_bootstrap.lme <- function(model, fn, B){
  # getVarCov.lme is the limiting factor...
  if (length(model$group) > 1) 
    stop("not implemented for multiple levels of nesting")
  
  # Match function
  fn <- match.fun(fn)
  
  # Extract fixed effects
  model.fixef <- nlme::fixef(model)
  
  ystar <- nlmeU::simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  ystar <- data.frame(ystar)
  
  t0 <- fn(model)
  
  # t.res <- list()
  # length(t.res) <- B
  #    t.res <- matrix(0, ncol = 2, nrow = B)
  #   for(i in 1:B){
  #     myin <- ystar[,i]
  #     model.update <- update(object = model, fixed = myin ~ .)
  #     t.res[i,] <- fn(model.update)
  #   }
  
  
  
  res <- lapply(ystar, function(y) {
    fit <- tryCatch(fn(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
#   for(i in 1:B){
#     #     myin <- ystar[,i]
#     #     model.update <- nlme:::update.lme(object = model, fixed = myin ~ .)
#     #     t.res[i,] <- fn(model.update)
#     fit <- tryCatch(fn(updated.model(model = model, new.y = ystar[,i])),  error = function(e) e)
#     if (inherits(fit, "error")) {
#       structure(rep(NA, length(t0)), fail.msgs = fit$message)
#     } else{
#       t.res[[i]] <- fit
#     }
    ### NOTE: Need to check if there was an error in try, and then use fn()
    # tmp.mod <- updated.model(model = model, new.y = ystar[,i])
    # t.res[[i]] <- fn(try.fit)
  # }
  tstar <- do.call('cbind', res)
  # tstar <- data.frame(tstar)
  
  #   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
  #   
  #   tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}



#' @rdname case_bootstrap
#' @export
case_bootstrap.lme <- function (model, fn, B, resample){
  
  data <- model$data
  # data$.id <- seq_len(nrow(data))
  clusters <- c(names(model$groups), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  t0 <- fn(model)
  
  rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(dat = data, cluster = clusters, resample = resample))))
  
  res <- lapply(rep.data, function(df) {
    fit <- tryCatch(fn(updated.model(model = model, new.data = df)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
  tstar <- do.call('cbind', res)
  
  rownames(tstar) <- names(t0)
  colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "case", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}


#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lme <- function (model, fn, B){
  fn <- match.fun(fn)
  
  t0 <- fn(model)
  ystar <- lapply(1:B, function(x) .resample.resids.lme(model))

#   return(.bootstrap.completion(model, ystar, B, fn))
  
  res <- lapply(ystar, function(y) {
    fit <- tryCatch(fn(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
  tstar <- do.call('cbind', res)

#   rownames(tstar) <- names(t0)
  colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "resid", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}


.resample.resids.lme <- function(model){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, level = 0) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- nlme::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  level.num <- ncol(model$groups)
  
  # Extract Zt (like lme4) design matrix
  re.form <- formula(model$modelStruct$reStr)
  Z <- lapply(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  
  if(level.num == 1) {
    bstar <- sample(model.ranef, replace = TRUE)
    
    Z <- as.data.frame(Z[[1]])
    Zlist <- lapply(Z, function(col) split(col, model$group))
    
    Zbstar <- lapply(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], bstar[,j], SIMPLIFY = FALSE) ))
    Zbstar.sum <- Reduce("+", Zbstar)
  } else{
    bstar <- lapply(model.ranef,
                    FUN = function(x) {
                      J <- nrow(x)
                      
                      # Sample of b*
                      bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                      bstar <- x[bstar.index,]
                      return(bstar)
                    })
    
    
    Z <- lapply(Z, function(zi) as.data.frame(zi))
    Z  <- Z [rev(names(Z))] # agree w/ order of model$group and bstar
    
    Zlist <- lapply(1:length(Z), function(i) lapply(Z[[i]], function(col) split(col, model$group[,i])))
    names(Zlist) <- names(Z)
    
    
    Zbstar <- lapply(1:length(Zlist), function(e) {
      z.e <- Zlist[[e]]
      b.e <- bstar[[e]]
      if(is.numeric(b.e)){
        unlist(lapply(1:length(e), function(j) unlist(mapply("*", z.e[[j]], b.e, SIMPLIFY = FALSE) )), recursive = FALSE)
      } else{
        unlist(lapply(1:length(e), function(j) unlist(mapply("*", z.e[[j]], b.e[,j], SIMPLIFY = FALSE) )), recursive = FALSE)
      }
    })
    
    Zbstar.sum <- Reduce("+", Zbstar)
    
    
#     Zlist <- lapply(Zlist, function(e) lapply(e, function(x) Matrix::bdiag(x)))
#     names(Zlist) <- names(Z)
#     
#     Zlist  <- unlist(Zlist, recursive = FALSE) 
#     
#     bstar <- bstar[rev(names(bstar))] # agree w/ order of Zlist
#     names(bstar) <- names(Zlist)
#     
#     Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Z)
#     Zbstar.sum <- Reduce("+", Zbstar)
    
  }
  
  
  # Resample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
  
}


#' @rdname reb_bootstrap
#' @export
reb_bootstrap.lme <- function (model, fn, B, reb_type = 0){
  
  if(ncol(model$groups) > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  if(reb_type != 2) fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.reb.lme(model = model, reb_type = reb_type)) )
  
 
  
  if(reb_type == 2){
    fe.0 <- nlme::fixef(model)
    vc.0 <- nlme::getVarCov(model)
    t0 <- c(beta = fe.0, sigma = c(diag(vc.0), model$sigma^2))
    tstar <- lapply(ystar, function(y) {
      fit <- tryCatch(updated.model(model = model, new.y = y),  
                      error = function(e) e)
      if (inherits(fit, "error")) {
        structure(list(poi = rep(NA, length(nlme::fixef(model))), varcomp = rep(NA, length(diag(vc.0)) + 1)), 
                  fail.msgs = fit$message)
      } else{
        vc <- nlme::getVarCov(fit)
        list(fixef = nlme::fixef(fit), varcomp = unname(c(diag(vc), fit$sigma^2)))
      }
    })
    
    vcs <- lapply(tstar, function(x) x$varcomp)
    Sb <- log( do.call("rbind", vcs) )
#     fes <- lapply(tstar, function(x) x$fixef)
    
    Mb <- matrix(rep(apply(Sb, 2, mean, na.rm = TRUE), times = B), nrow = B, byrow = TRUE)
    CovSb <- cov(na.omit(Sb))
    SdSb <- sqrt(diag(CovSb))
    
    Db <- matrix(rep(SdSb, times = B), nrow = B, byrow = TRUE)
    
    EW <- eigen(solve(CovSb), symmetric = T)
    Whalf <- EW$vectors %*% diag(sqrt(EW$values))
    
    Sbmod <- (Sb - Mb) %*% Whalf
    Sbmod <- Sbmod * Db # elementwise not a type 
    Lb <- exp(Mb + Sbmod)
    
    tstar <- lapply(tstar, unlist)
  } else{
    t0 <- fn(model)
#     Lb <- NULL
    tstar <- lapply(ystar, function(y) {
      fit <- tryCatch(fn(updated.model(model = model, new.y = y)),  
                      error = function(e) e)
      if (inherits(fit, "error")) {
        structure(rep(NA, length(t0)), fail.msgs = fit$message)
      } else{
        fit
      }
    }
    )
  }
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
#   rownames(tstar) <- names(fn(model))
  
  if(reb_type == 2) {
    idx <- 1:length(fe.0)
    fe.star <- tstar[idx,] 
    fe.adj <- sweep(fe.star, MARGIN = 1, STATS = fe.0 - rowMeans(fe.star, na.rm = TRUE), FUN = "+")
      
    vc.star <- tstar[-idx,] 
    vc.adj <- sweep(vc.star, MARGIN = 1, STATS = t0[-idx] / rowMeans(vc.star, na.rm = TRUE), FUN = "*")
    
    tstar <- rbind(fe.adj, vc.adj)
    
    fn <- function(.) {
      c(beta = nlme::fixef(.), sigma = c(diag(nlme::getVarCov(.)), .$sigma^2))
    }
  }

  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = paste("reb", reb_type, sep = ""), call = match.call()),
                   class = "boot")
  
  return(RES)
  
}


#' REB resampling procedures
#' @importFrom RLRsim extract.lmeDesign
.resample.reb.lme <- function(model, reb_type){
  
  dsgn <- RLRsim::extract.lmeDesign(model)
  
  # extract marginal residuals
  model.mresid <- dsgn$y - predict(model, level = 0)
  
  # Extract Z design matrix
  Z <- Matrix::Matrix(dsgn$Z)
  
  # level 2 resid
  u <- solve(t(Z) %*% Z) %*% t(Z) %*% model.mresid # a single vector
  
  # level 1 resid
  e <- model.mresid - Z %*% u
  
  # The current way u is organized is inspired by the 
  # ranef.merMod function in lme4.
  # TODO: think about 3+ level models...
  #   ans <- model@pp$b(1)
  levs <- lapply(fl <- model$groups, levels)
  asgn <- seq_along(fl)
  re <- model$coefficients$random
  cnms <- lapply(re, colnames)
  nc <- vapply(cnms, length, 1L)
  nb <- nc * (nl <- vapply(levs, length, 1L))
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(u, nbseq)
  for (i in seq_along(u))
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  names(u) <- names(cnms)
  
  if(reb_type == 1){
    # Calculations
    Uhat <- lapply(u, function(x){
      S <- (t(x) %*% x) / nrow(x)
      R <- nlme::getVarCov(model)
      Ls <- chol(S, pivot = TRUE)
      Lr <- chol(R, pivot = TRUE)
      A <- t(Lr %*% solve(Ls))
      
      Uhat <- x%*%A
      
      # center
      Uhat <- data.frame(scale(Uhat, scale = FALSE))
      
      return(Uhat)
    })
    
    sigma <- model$sigma
    estar <- sigma * e %*% ((t(e) %*% e) / length(e))^(-1/2)
    estar <- data.frame(scale(estar, scale = FALSE))
    
  } else{
    Uhat <- u
    estar <- e
  }
  
  Xbeta <- predict(model, level = 0)
  
  # resample uhats
  ustar <- lapply(Uhat,
                  FUN = function(x) {
                    J <- nrow(x)
                    x <- as.data.frame(x)
                    # Sample of b*
                    ustar <- sample(x, replace = TRUE)
                    return(ustar)
                  })
   
  ## TODO: fix this issue with the levels... Need to resample from here...
  
  # Extract Z design matrix separated by variance
  re.form <- formula(model$modelStruct$reStr)
  Z <- lapply(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  Z <- as.data.frame(Z[[1]])
  Zlist <- lapply(Z, function(col) split(col, model$group))

  ustar <- ustar[[1]] # since only working with 2-levels models now
  
  # Get Zb*
  Zbstar <- lapply(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], ustar[,j], SIMPLIFY = FALSE) ))
  Zbstar.sum <- Reduce("+", Zbstar)

  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}


#' @rdname cgr_bootstrap
#' @export
cgr_bootstrap.lme <- function (model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.cgr.lme(model = model)) )
  
  t0 <- fn(model)
  
  
  res <- lapply(ystar, function(y) {
    fit <- tryCatch(fn(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
  tstar <- do.call("cbind", res) # Can these be nested?
  colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "CGR", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' CGR resampling procedures
.resample.cgr.lme <- function(model){
  level.num <- ncol(model$groups)
  
  # Extract random effects
  model.ranef <- nlme::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Higher levels
  if(level.num == 1) {
    model.ranef <- list(model.ranef)
    names(model.ranef) <- colnames(model$groups)
  }
  
  re.struct <- model$modelStruct$reStruct
  sigma <- model$sigma
  
  Uhat.list <- lapply(seq_along(model.ranef),
                      FUN = function(i) {
                        u <- scale(model.ranef[[i]], scale = FALSE)
                        S <- (t(u) %*% u) / length(u)
                        
                        re.name <- names(model.ranef)[i]
                        R <- sigma^2 * as.matrix(re.struct[[re.name]])
                        
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
  ehat <- sigma * e * ((t(e) %*% e) / length(e))^(-1/2)
  
  # Resample Uhat
  ustar <- lapply(Uhat.list,
                  FUN = function(df) {
                    index <- sample(x = seq_len(nrow(df)), size = nrow(df), replace = TRUE)
                    return(as.data.frame(df[index,]))
                  })
  
  # Extract Z design matrix
  # Extract Zt (like lme4) design matrix
  re.form <- formula(re.struct)
  Z <- lapply(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  
  if(level.num == 1) {
    ustar <- ustar[[1]]
    
    Z <- as.data.frame(Z[[1]])
    Zlist <- lapply(Z, function(col) split(col, model$group))
    
    Zbstar <- lapply(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], ustar[,j], SIMPLIFY = FALSE) ))
  } else{
    
    Z <- lapply(Z, function(zi) as.data.frame(zi))
    Z  <- Z [rev(names(Z))] # agree w/ order of model$group and bstar
    
    Zlist <- lapply(1:length(Z), function(i) lapply(Z[[i]], function(col) split(col, model$group[,i])))
    names(Zlist) <- names(Z)
    
    
    Zbstar <- lapply(1:length(Zlist), function(e) {
      z.e <- Zlist[[e]]
      u.e <- ustar[[e]]
      if(is.numeric(u.e)){
        unlist(lapply(1:length(e), function(j) unlist(mapply("*", z.e[[j]], u.e, SIMPLIFY = FALSE) )), recursive = FALSE)
      } else{
        unlist(lapply(1:length(e), function(j) unlist(mapply("*", z.e[[j]], u.e[,j], SIMPLIFY = FALSE) )), recursive = FALSE)
      }
    })
  }
  
  Zbstar.sum <- Reduce("+", Zbstar)  
  
  # Extract fixed part of the model
  Xbeta <- predict(model, level = 0)
  
  # Get e*
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  # Combine
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}