#' @rdname bootstrap
#' @export
#' @importFrom stats as.formula cov formula model.matrix na.exclude na.omit predict resid simulate
bootstrap.lmerMod <- function (model, fn, type, B, resample, reb_type){
  switch(type,
         parametric = parametric_bootstrap.lmerMod(model, fn, B),
         residual = resid_bootstrap.lmerMod(model, fn, B),
         case = case_bootstrap.lmerMod(model, fn, B, resample),
         cgr = cgr_bootstrap.lmerMod(model, fn, B),
         reb = reb_bootstrap.lmerMod(model, fn, B, reb_type = 0))
  # TODO: need to be able to save results
}


#' @rdname parametric_bootstrap
#' @export
parametric_bootstrap.lmerMod <- function(model, fn, B){
  fn <- match.fun(fn)

  model.fixef <- lme4::fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)

  return(.bootstrap.completion(model, ystar, B, fn))

  # TODO: once we have things working, think about parallelization.
  #       using an llply statement would make this easy with the .parallel
  #       parameter, but it might be slower than using mclapply, which is
  #       found in the parallel package.
}


#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lmerMod <- function (model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- lapply(1:B, function(x) .resample.resids(model))

  RES <- .bootstrap.completion(model, ystar, B, fn)
  RES$sim <- "resid"
  return(RES)
}


#' @rdname case_bootstrap
#' @export
case_bootstrap.lmerMod <- function (model, fn, B, resample){
  
  data <- model@frame
  # data$.id <- seq_len(nrow(data))
  clusters <- c(rev(names(lme4::getME(model, "flist"))), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  rep.data <- lapply(integer(B), function(x) .cases.resamp(dat = data, cluster = clusters, resample = resample))
  
  # Plugin to .cases.completion due to small changes
  RES <- .cases.completion(model, rep.data, B, fn)
  return(RES)
}


# # Using recursion allows for a very general function...
# # How can I speed this up?
# .cases.resamp <- function(dat, cluster, resample) {
#   # exit early for trivial data
#   if(nrow(dat) == 1 || all(resample==FALSE))
#     return(dat)
#   
#   # sample the clustering factor
#   cls <- sample(unique(dat[[cluster[1]]]), resample=resample[1])
#   
#   # subset on the sampled clustering factors
#   sub <- lapply(cls, function(b) dat[dat[[cluster[1]]]==b,])
#   
#   # sample lower levels of hierarchy (if any)
#   if(length(cluster) > 1)
#     sub <- lapply(sub, .cases.resamp, cluster=cluster[-1], resample=resample[-1])
#   
#   # join and return samples
#   do.call(rbind, sub)
# }
# 
.cases.resamp <- function(dat, cluster, resample) {
  # exit early for trivial data
  if(nrow(dat) == 1 || all(resample==FALSE))
    return(dat)
  
  # ver <- as.numeric_version(packageVersion("dplyr"))
  res <- dat
  
  for(i in 1:length(cluster)) {
    
    if(i==1 & resample[i]) {
      dots <- cluster[1]
      grouped <- dplyr::group_by_(res, dots)
      g_rows <- dplyr::group_rows(grouped)
      # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- cluster[-i]
        grouped <- dplyr::group_by_(res, .dots = dots)
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- cluster[i]
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- plyr::ldply(res, function(df) {
            grouped <- dplyr::group_by_(df, .dots = dots)
            g_rows <- dplyr::group_rows(grouped)
            # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
            cls <- sample(seq_along(g_rows), replace = resample[i])
            idx <- unlist(g_rows[cls], recursive = FALSE)
            grouped[idx, ]
          }, .id = NULL)
        }
      }
    }
    
    
  }
  return(res)
}


# .cases.resamp <- function (model, extra_step){
#   # Draw sample of size J from level-2 units
#   model.split <- split(x=model@frame, f=model@flist)
#   model.split.samp <- sample(x=model.split, size = length(model.split), resample = TRUE)
#   # For each sample, draw a sample of the cases from the level-2 unit
#   
#   if(extra_step == TRUE){
#     model.resamp <- lapply(model.split.samp,
#                            FUN = function(x) {
#                              J <- nrow(x)
#                              
#                              # Sample of level-2 rows
#                              model.sub.index <- sample(x = seq_len(J), size = J, resample = TRUE)
#                              resampled <- x[model.sub.index,]
#                              return(resampled)
#                            })
#     model.comb <- do.call('rbind', model.resamp)
#   } else{ # else statement needs to be located here
#     model.comb <- do.call('rbind', model.split.samp)
#   }
# }

.cases.completion <- function(model, data, B, fn){
  t0 <- fn(model)
  
  # Refit the model and apply 'fn' to it using lapply
  form <- model@call$formula
  reml <- lme4::isREML(model)
  tstar <- lapply(data, function(x) {
    fn(lme4::lmer(formula = form, data = x, REML = reml)) 
  })
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "case", call = match.call()),
                   class = "boot")
  
  return(RES)
}


#' @rdname cgr_bootstrap
#' @export
cgr_bootstrap.lmerMod <- function (model, fn, B){
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.cgr(model = model)) )
  
  RES <- .bootstrap.completion(model, ystar, B, fn)
  RES$sim <- "cgr"
  return(RES)
  
}


#' @rdname reb_bootstrap
#' @export
reb_bootstrap.lmerMod <- function (model, fn, B, reb_type = 0){
  
  if(lme4::getME(object = model, name = "n_rfacs") > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  fn <- match.fun(fn)
  
  ystar <- as.data.frame( replicate(n = B, .resample.reb(model = model, reb_type = reb_type)) )
  
  if(reb_type != 2) t0 <- fn(model)
  
  # Refit the model and apply 'fn' to it using lapply
#   tstar <- lapply(ystar[1,], function(x) {
#     fn(refit(object = model, newresp = x))
#   })

  if(reb_type == 2){
    fe.0 <- lme4::fixef(model)
    vc.0 <- bdiag(lme4::VarCorr(model))
    t0 <- c(beta = fe.0, sigma = c(diag(vc.0), lme4::getME(model, "sigma")^2))
    tstar <- lapply(ystar, function(x) {
      m <- lme4::refit(object = model, newresp = x)
      vc <- as.data.frame(lme4::VarCorr(m))
      list(fixef = lme4::fixef(m), varcomp = vc$vcov[is.na(vc$var2)])
    })
    
    vcs <- lapply(tstar, function(x) x$varcomp)
    Sb <- log( do.call("rbind", vcs) )
#     fes <- lapply(tstar, function(x) x$fixef)
    
    Mb <- matrix(rep(apply(Sb, 2, mean), times = B), nrow = B, byrow = TRUE)
    CovSb <- cov(Sb)
    SdSb <- sqrt(diag(CovSb))
    
    Db <- matrix(rep(SdSb, times = B), nrow = B, byrow = TRUE)
    
    EW <- eigen(solve(CovSb), symmetric = T)
    Whalf <- EW$vectors %*% diag(sqrt(EW$values))

    Sbmod <- (Sb - Mb) %*% Whalf
    Sbmod <- Sbmod * Db # elementwise not a type 
    Lb <- exp(Mb + Sbmod)
    
    tstar <- lapply(tstar, unlist)
  } else{
    Lb <- NULL
    tstar <- lapply(ystar, function(x) {
      fn(lme4::refit(object = model, newresp = x))
    })
  }

  tstar <- do.call("cbind", tstar) # Can these be nested?
#   rownames(tstar) <- names(fn(model))


  if(reb_type == 2) {
    idx <- 1:length(fe.0)
    fe.star <- tstar[idx,] 
    fe.adj <- sweep(fe.star, MARGIN = 1, STATS = fe.0 - rowMeans(fe.star, na.rm = TRUE), FUN = "+")
    
    vc.star <- tstar[-idx,] 
    vc.adj <- sweep(vc.star, MARGIN = 1, STATS = t0[-idx] / rowMeans(vc.star, na.rm = TRUE), FUN = "*")
    
    tstar <- rbind(fe.adj, vc.adj)
    
    fn <- function(.) {
      c(beta = lme4::fixef(.), sigma =c(diag(bdiag(lme4::VarCorr(.))), lme4::getME(., "sigma")^2))
    }
  }

    
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "reb", call = match.call(), reb2 = Lb),
                   class = "boot")
  
  return(RES)
}


#####################
# Utility Functions #
#####################

#' @title Zbstar combine
#'
#' @description
#' Combine \code{bstar} and \code{zstar} to create {Zbstar}.
#'
#' @details
#' This function combines \code{bstar} and \code{zstar} to create {Zbstar} using an lapply statement.
#'
#' @param bstar A list of matrices bstar
#' @param zstar A list of matrices zstar
#'
#' @return matrix
#' @keywords internal
#' @noRd
.Zbstar.combine <- function(bstar, zstar){
  lapply(1:length(bstar), function(i){
    Matrix::t(zstar[[i]]) %*% bstar[[i]]
  })
}


#' @title Bootstrap Completion
#'
#' @description
#' Finishes the bootstrap process and makes the output readable.
#'
#' @details
#' This function is given \code{model, ystar, B, fn} and uses them to complete
#' the bootstrap process. They are then structured into a list for output and returned.
#'
#' @param ystar The ystar being passed in
#' @inheritParams bootstrap
#'
#' @return list
#' @keywords internal
#' @noRd
.bootstrap.completion <- function(model, ystar, B, fn){
  t0 <- fn(model)

  # Refit the model and apply 'fn' to it using lapply
  tstar <- lapply(ystar, function(x) {
    fn(lme4::refit(object = model, newresp = x))
  })

  nsim <- length(tstar)
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  if ((nfail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", nfail, "/", nsim, ")")
    fail.msgs <- vapply(tstar[bad.runs], FUN=attr, FUN.VALUE = character(1),
                        "fail.msgs")
  } else fail.msgs <- NULL

  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  attr(RES,"bootFail") <- nfail
  attr(RES,"boot.fail.msgs") <- fail.msgs
  attr(RES,"boot_type") <- "boot"
  return(RES)
}

#' CGR resampling procedures
#' @keywords internal
#' @noRd
.resample.cgr <- function(model){
  model.ranef <- lme4::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Higher levels
  Uhat.list <- lapply(seq_along(model.ranef),
                      FUN = function(i) {
                        u <- scale(model.ranef[[i]], scale = FALSE)
                        S <- (t(u) %*% u) / length(u)
                        
                        re.name <- names(model.ranef)[i]
                        R <- bdiag(lme4::VarCorr(model)[[names(model.ranef)[i]]])
                        
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
  sigma <- lme4::getME(model, "sigma")
  ehat <- sigma * e * as.numeric((t(e)%*%e) / length(e))^(-1/2)
  
  # Extract Z design matrix
  Z <- lme4::getME(object = model, name = "Ztlist")
  
  
  Xbeta <- predict(model, re.form = NA)
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  # Resample Uhat
  ustar <- lapply(Uhat.list,
                  FUN = function(df) {
                    index <- sample(x = seq_len(nrow(df)), size = nrow(df), replace = TRUE)
                    return(df[index,])
                  })
  
  # Structure u*
  if(level.num == 1){
    if(is.data.frame(ustar[[1]])){
      ustar <- lapply(ustar, FUN = function(x) as.list(x))[[1]] 
    }
    names(ustar) <- names(Z)
  } else {
    ustar <- lapply(ustar, FUN = function(x) as.data.frame(x))
    ustar <- do.call(c, ustar)
    names(ustar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Get e*
  estar <- sample(x = ehat, size = length(ehat), replace = TRUE)
  
  # Combine
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}

#' Resampling residuals from mixed models
#' @keywords internal
#' @noRd
.resample.resids <- function(model){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- lme4::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Extract Z design matrix
  Z <- lme4::getME(object = model, name = "Ztlist")
  
  bstar <- lapply(model.ranef,
                  FUN = function(x) {
                    J <- nrow(x)
                    
                    # Sample of b*
                    bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                    bstar <- x[bstar.index,]
                    return(bstar)
                  })
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    if(!is.numeric(bstar[[1]])) bstar <- lapply(bstar, FUN = function(x) as.list(x))[[1]]
    names(bstar) <- names(Z)
  } else {
    bstar <- lapply(bstar, FUN = function(x) as.data.frame(x))
    bstar <- do.call(c, bstar)
    names(bstar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  
  # Resample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
  
}

#' REB resampling procedures
#' 
#' @param reb_type Specifies the inclusion of REB/1
#' @inheritParams bootstrap
#' @import Matrix
#' @keywords internal
#' @noRd
.resample.reb <- function(model, reb_type){
  # extract marginal residuals
  model.mresid <- lme4::getME(model, "y") - predict(model, re.form = NA)
  
  # Extract Z design matrix
  Z <- lme4::getME(object = model, name = "Z")
  
  # level 2 resid
  u <- solve(t(Z) %*% Z) %*% t(Z) %*% model.mresid # a single vector
  
  # level 1 resid
  e <- model.mresid - Z %*% u
  
  # The current way u is organized is inspired by the 
  # ranef.merMod function in lme4.
  # TODO: think about 3+ level models...
  #   ans <- model@pp$b(1)
  levs <- lapply(fl <- model@flist, levels)
  asgn <- attr(fl, "assign")
  cnms <- model@cnms
  nc <- vapply(cnms, length, 1L)
  nb <- nc * (nl <- vapply(levs, length, 1L)[asgn])
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(u, nbseq)
  for (i in seq_along(u))
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  names(u) <- names(cnms)
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  if(reb_type == 1){
    if(level.num > 1) stop("reb_type = 1 is not yet implemented for higher order models")
    # Calculations
    Uhat <- lapply(u, function(x){
      S <- (t(x) %*% x) / nrow(x)
      R <- bdiag(lme4::VarCorr(model))
      Ls <- chol(S, pivot = TRUE)
      Lr <- chol(R, pivot = TRUE)
      A <- t(Lr %*% solve(Ls))
      
      Uhat <- x%*%A
      
      # center
      Uhat <- data.frame(scale(Uhat, scale = FALSE))
      
      return(Uhat)
    })
    
    sigma <- lme4::getME(model, "sigma")
    estar <- sigma * e %*% ((t(e) %*% e) / length(e))^(-1/2)
    estar <- data.frame(scale(estar, scale = FALSE))
    
  } else{
    Uhat <- u
    estar <- e
  }
  
  Xbeta <- predict(model, re.form = NA)
  
  # resample uhats
  ustar <- lapply(Uhat,
                  FUN = function(x) {
                    J <- nrow(x)
                    
                    # Sample of b*
                    ustar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                    ustar <- data.frame(x[ustar.index,])
                    return(ustar)
                  })
  
#   Uhat <- as.data.frame(as.matrix(Uhat))
#   Uhat.list <- list(Uhat)
  
  ## TODO: fix this issue with the levels... Need to resample from here...
  
  # Extract Z design matrix separated by variance
  Ztlist <- lme4::getME(object = model, name = "Ztlist")
  
  if(level.num == 1){
    ustar <- lapply(ustar, FUN = function(x) as.list(x))[[1]]
  } else {
    ustar <- lapply(ustar, FUN = function(x) as.data.frame(x))
    ustar <- do.call(c, ustar)
  }
  
  names(ustar) <- names(Ztlist) 
  ustar.df <- as.data.frame(ustar)
  
#   if(level.num == 1){
#     Uhat.list <- lapply(Uhat.list, FUN = function(x) as.list(x))[[1]]
#     names(Uhat.list) <- names(Ztlist)
#   } else {
#     Uhat.list <- sapply(Uhat.list, FUN = function(x) as.list(x))
#   }
  
#   # Resample Uhat
#   ustar <- sample(x = Uhat.list[[1]], size = length(Uhat.list[[1]]), replace = TRUE)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar.df, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = model.mresid, size = length(model.mresid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
}
