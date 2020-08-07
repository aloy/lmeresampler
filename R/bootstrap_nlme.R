#' @rdname bootstrap
#' @export
bootstrap.lme <- function(model, .f, type, B, resample, reb_type, linked){
  switch(type,
         parametric = parametric_bootstrap.lme(model, .f, B, type = type),
         residual = resid_bootstrap.lme(model, .f, B, type = type, linked = FALSE),
         case = case_bootstrap.lme(model, .f, B, resample, type = type),
         cgr = cgr_bootstrap.lme(model, .f, B, type = type),
         reb = reb_bootstrap.lme(model, .f, B, reb_type = 0))
}


#' @rdname parametric_bootstrap
#' @export
#' @importFrom nlmeU simulateY
parametric_bootstrap.lme <- function(model, .f, B, type){
  # getVarCov.lme is the limiting factor...
  
  # if(ncol(model$groups) > 1){
  #   stop("The REB bootstrap has not been adapted for 3+ level models.")
  # }
  
  # Match function
  .f <- match.fun(.f)
  
  # Extract fixed effects
  model.fixef <- nlme::fixef(model)
  
  ystar <- nlmeU::simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  ystar <- data.frame(ystar)
  
  t0 <- .f(model)
  
  # t.res <- list()
  # length(t.res) <- B
  #    t.res <- matrix(0, ncol = 2, nrow = B)
  #   for(i in 1:B){
  #     myin <- ystar[,i]
  #     model.update <- update(object = model, fixed = myin ~ .)
  #     t.res[i,] <- .f(model.update)
  #   }
  
  tstar <- purrr::map(ystar, function(y) {
    fit <- tryCatch(.f(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if(inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  })
  
  #   for(i in 1:B){
  #     #     myin <- ystar[,i]
  #     #     model.update <- nlme:::update.lme(object = model, fixed = myin ~ .)
  #     #     t.res[i,] <- .f(model.update)
  #     fit <- tryCatch(.f(updated.model(model = model, new.y = ystar[,i])),  error = function(e) e)
  #     if (inherits(fit, "error")) {
  #       structure(rep(NA, length(t0)), fail.msgs = fit$message)
  #     } else{
  #       t.res[[i]] <- fit
  #     }
  ### NOTE: Need to check if there was an error in try, and then use .f()
  # tmp.mod <- updated.model(model = model, new.y = ystar[,i])
  # t.res[[i]] <- .f(try.fit)
  # }
  tstar <- do.call('cbind', tstar)
  # tstar <- data.frame(tstar)
  
  #   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
  #   
  #   tstar <- do.call("cbind", tstar) # Can these be nested?
  row.names(tstar) <- names(t0)
  # colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  if((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- purrr::map_chr(tstar[bad.runs], .f = function(x){
      attr(x)})
  } else fail.msgs <- NULL 
  
  # prep for stats df
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model$data,
                        seed = .Random.seed, type = type, call = match.call()), 
                   class = "lmeresamp")
  
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' @rdname case_bootstrap
#' @export
case_bootstrap.lme <- function(model, .f, B, resample, type){
  
  data <- model$data
  # data$.id <- seq_len(nrow(data))
  clusters <- c(names(model$groups), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  # rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(dat = data, cluster = clusters, resample = resample))))
  tstar <- purrr::map(integer(B), function(x) .cases.resamp(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  # res <- purrr::map(rep.data, function(df) {
  #   fit <- tryCatch(.f(updated.model(model = model, new.data = df)),  
  #                   error = function(e) e)
  #   if (inherits(fit, "error")) {
  #     structure(rep(NA, length(t0)), fail.msgs = fit$message)
  #   } else{
  #     fit
  #   }
  # })
  
  t0 <- .f(model)
  tstar <- do.call("cbind", tstar)
  row.names(tstar) <- names(t0)
  # colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  if((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- purrr::map_chr(tstar[bad.runs], .f = attr,  FUN.VALUE = character(1),
                                "fail.msgs")
  } else fail.msgs <- NULL
  
  # prep for stats df
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model$data,
                        seed = .Random.seed, type = type, call = match.call()),
                   class = "lmeresamp")
  
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lme <- function(model, .f, B, type, linked = FALSE){
  .f <- match.fun(.f)
  
  t0 <- .f(model)
  tstar <- purrr::map(1:B, function(x) .resample.resids.lme(model, .f, linked = linked))
  
  tstar <- do.call('cbind', tstar)
  
  #   rownames(tstar) <- names(t0)
  # colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- purrr::map_chr(res[bad.runs], .f = attr,  FUN.VALUE = character(1),
                                "fail.msgs")
  } else fail.msgs <- NULL
  
  # prep for stats df
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model$data,
                        seed = .Random.seed, type = type, call = match.call()),
                   class = "lmeresamp")
  
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' @keywords internal
#' @noRd
.resample.resids.lme <- function(model, .f, linked = linked){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, level = 0) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- nlme::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  level.num <- ncol(model$groups)
  
  # Extract Zt (like lme4) design matrix
  re.form <- formula(model$modelStruct$reStr)
  Z <- purrr::map(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  
  if(level.num == 1) {
    bstar <- sample(model.ranef, replace = TRUE)
    
    Z <- as.data.frame(Z[[1]])
    Zlist <- purrr::map(Z, function(col) split(col, model$group))
    
    # purrr::map2, check output and see what each apply does
    Zbstar <- purrr::map(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], bstar[,j], SIMPLIFY = FALSE)))
    Zbstar.sum <- Reduce("+", Zbstar)
  } else{
    bstar <- purrr::map(model.ranef,
                        .f = function(x) {
                          J <- nrow(x)
                          
                          # Sample of b*
                          bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                          bstar <- x[bstar.index,]
                          return(bstar)
                        })
    
    Z <- purrr::map(Z, function(zi) as.data.frame(zi))
    Z  <- Z [rev(names(Z))] # agree w/ order of model$group and bstar
    
    Zlist <- purrr::map(1:length(Z), function(i) purrr::map(Z[[i]], function(col) split(col, model$group[,i])))
    names(Zlist) <- names(Z)
    
    Zbstar <- purrr::map(1:length(Zlist), function(e) {
      z.e <- Zlist[[e]]
      b.e <- bstar[[e]]
      if(is.numeric(b.e)){
        unlist(purrr::map(1:length(e), function(j) unlist(mapply("*", z.e[[j]], b.e, SIMPLIFY = FALSE) )), recursive = FALSE)
      } else{
        unlist(purrr::map(1:length(e), function(j) unlist(mapply("*", z.e[[j]], b.e[,j], SIMPLIFY = FALSE) )), recursive = FALSE)
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
  
  if(linked == FALSE){
    # Resample residuals
    estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
    
    # Combine function
    y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
    
    # .f(model) is t0
    tstar <- tryCatch(.f(updated.model(model = model, new.y = y.star)),  
                      error = function(e) e)
    if(inherits(tstar, "error")) {
      structure(rep(NA, length(.f(model))), fail.msgs = tstar$message)
    } else{
      tstar
    }
  } else {
    #linked
    model.mresid <- nlme::getResponse(model) - predict(model, re.form = ~0)
    model.mresid.cent <- scale(model.mresid, scale = FALSE)
    
    # Resample residuals
    mresid.star <- sample(x = model.mresid.cent, size = length(model.mresid.cent), replace = TRUE)
    
    # Combine function
    y.star <- as.numeric(Xbeta + mresid.star)
    
    # .f(model) is t0
    # Refit
      tstar <- tryCatch(.f(updated.model(model = model, new.y = y.star)),  
                      error = function(e) e)
      if(inherits(tstar, "error")) {
        structure(rep(NA, length(.f(model))), fail.msgs = tstar$message)
      } else{
        tstar
      }
  }
  return(tstar)
}


#' @rdname reb_bootstrap
#' @inheritParams bootstrap
#' @export
reb_bootstrap.lme <- function(model, .f, B, reb_type = 0){
  
  if(ncol(model$groups) > 1){
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  if(reb_type != 2) .f <- match.fun(.f)
  
  ystar <- as.data.frame(replicate(n = B, .resample.reb.lme(model = model, .f, reb_type = reb_type)))
  
  if(reb_type == 2){
    fe.0 <- nlme::fixef(model)
    vc.0 <- nlme::getVarCov(model)
    t0 <- c(beta = fe.0, sigma = c(diag(vc.0), model$sigma^2))
    tstar <- purrr::map(ystar, function(y) {
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
    
    vcs <- purrr::map(tstar, function(x) x$varcomp)
    Sb <- log(do.call("rbind", vcs))
    #     fes <- lapply(tstar, function(x) x$fixef)
    
    Mb <- matrix(rep(colMeans(Sb, na.rm = TRUE), times = B), nrow = B, byrow = TRUE)
    CovSb <- cov(na.omit(Sb))
    SdSb <- sqrt(diag(CovSb))
    
    Db <- matrix(rep(SdSb, times = B), nrow = B, byrow = TRUE)
    
    EW <- eigen(solve(CovSb), symmetric = T)
    Whalf <- EW$vectors %*% diag(sqrt(EW$values))
    
    Sbmod <- (Sb - Mb) %*% Whalf
    Sbmod <- Sbmod * Db # elementwise not a type 
    Lb <- exp(Mb + Sbmod)
    
    tstar <- purrr::map(tstar, unlist)
  } else{
    t0 <- .f(model)
    #     Lb <- NULL
    tstar <- purrr::map(ystar, function(y) {
      fit <- tryCatch(.f(updated.model(model = model, new.y = y)),  
                      error = function(e) e)
      if (inherits(fit, "error")) {
        structure(rep(NA, length(t0)), fail.msgs = fit$message)
      } else{
        fit
      }
    })
  }
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  # colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
  #   rownames(tstar) <- names(.f(model))
  
  if(reb_type == 2) {
    idx <- 1:length(fe.0)
    fe.star <- tstar[idx,] 
    fe.adj <- sweep(fe.star, MARGIN = 1, STATS = fe.0 - rowMeans(fe.star, na.rm = TRUE), FUN = "+")
    
    vc.star <- tstar[-idx,] 
    vc.adj <- sweep(vc.star, MARGIN = 1, STATS = t0[-idx] / rowMeans(vc.star, na.rm = TRUE), FUN = "*")
    
    tstar <- rbind(fe.adj, vc.adj)
    
    .f <- function(.) {
      c(beta = nlme::fixef(.), sigma = c(diag(nlme::getVarCov(.)), .$sigma^2))
    }
  }
  
  # prep for stats df
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model$data,
                        seed = .Random.seed, type = paste("reb", reb_type, sep = ""), call = match.call()),
                   class = "lmeresamp")
  return(RES)
}


#' REB resampling procedures
#' @importFrom RLRsim extract.lmeDesign
#' @keywords internal
#' @noRd
.resample.reb.lme <- function(model, .f, reb_type){
  
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
  levs <- purrr::map(fl <- model$groups, levels)
  asgn <- seq_along(fl)
  re <- model$coefficients$random
  cnms <- purrr::map(re, colnames)
  nc <- purrr::map_int(cnms, length) #map_int()
  nb <- nc * (nl <- purrr::map_int(levs, length))
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(u, nbseq)
  for (i in seq_along(u))
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  names(u) <- names(cnms)
  
  if(reb_type == 1){
    # Calculations
    Uhat <- purrr::map(u, function(x){
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
  ustar <- purrr::map(Uhat,
                      .f = function(x) {
                        J <- nrow(x)
                        x <- as.data.frame(x)
                        # Sample of b*
                        ustar <- sample(x, replace = TRUE)
                        return(ustar)
                      })
  
  ## TODO: fix this issue with the levels... Need to resample from here...
  
  # Extract Z design matrix separated by variance
  re.form <- formula(model$modelStruct$reStr)
  Z <- purrr::map(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  Z <- as.data.frame(Z[[1]])
  Zlist <- purrr::map(Z, function(col) split(col, model$group))
  
  ustar <- ustar[[1]] # since only working with 2-levels models now
  
  # Get Zb*
  Zbstar <- purrr::map(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], ustar[,j], SIMPLIFY = FALSE)))
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}


#' @rdname cgr_bootstrap
#' @inheritParams bootstrap
#' @export
cgr_bootstrap.lme <- function(model, .f, B, type = type){
  .f <- match.fun(.f)
  
  tstar <- as.data.frame(replicate(n = B, .resample.cgr.lme(model = model, .f)))
  
  t0 <- .f(model)
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  # colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- purrr::map_chr(res[bad.runs], .f = attr,  FUN.VALUE = character(1),
                                "fail.msgs")
  } else fail.msgs <- NULL
  
  # prep for stats df
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model$data,
                        seed = .Random.seed, type = type, call = match.call()),
                   class = "lmeresamp")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' CGR resampling procedures
#' @keywords internal
#' @noRd
.resample.cgr.lme <- function(model, .f){
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
  
  Uhat.list <- purrr::map(seq_along(model.ranef),
                          .f = function(i) {
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
  ehat <- sigma * e * as.numeric((t(e) %*% e) / length(e))^(-1/2)
  
  # Resample Uhat
  ustar <- purrr::map(Uhat.list,
                      .f = function(df) {
                        index <- sample(x = seq_len(nrow(df)), size = nrow(df), replace = TRUE)
                        return(as.data.frame(df[index,]))
                      })
  
  # Extract Z design matrix
  # Extract Zt (like lme4) design matrix
  re.form <- formula(re.struct)
  Z <- purrr::map(1:length(re.form), function(i) model.matrix(formula(model$modelStruct$reStr)[[i]], data=model$data))
  names(Z) <- names(re.form)
  
  if(level.num == 1) {
    ustar <- ustar[[1]]
    
    Z <- as.data.frame(Z[[1]])
    Zlist <- purrr::map(Z, function(col) split(col, model$group))
    
    Zbstar <- purrr::map(1:length(Zlist), function(j) unlist(mapply("*", Zlist[[j]], ustar[,j], SIMPLIFY = FALSE) ))
  } else{
    
    Z <- purrr::map(Z, function(zi) as.data.frame(zi))
    Z  <- Z [rev(names(Z))] # agree w/ order of model$group and bstar
    
    Zlist <- purrr::map(1:length(Z), function(i) purrr:map(Z[[i]], function(col) split(col, model$group[,i])))
    names(Zlist) <- names(Z)
    
    
    Zbstar <- purrr::map(1:length(Zlist), function(e) {
      z.e <- Zlist[[e]]
      u.e <- ustar[[e]]
      if(is.numeric(u.e)){
        unlist(purrr::map(1:length(e), function(j) unlist(mapply("*", z.e[[j]], u.e, SIMPLIFY = FALSE))), recursive = FALSE)
      } else{
        unlist(purrr::map(1:length(e), function(j) unlist(mapply("*", z.e[[j]], u.e[,j], SIMPLIFY = FALSE))), recursive = FALSE)
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
  
  # Refit
  tstar <- tryCatch(.f(updated.model(model = model, new.y = y.star)),  
                    error = function(e) e)
  if (inherits(tstar, "error")) {
    structure(rep(NA, length(.f(model))), fail.msgs = tstar$message)
  } else{
    tstar
  }
  return(tstar)
}

