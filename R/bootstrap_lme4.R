#' @rdname bootstrap
#' @export
#' @importFrom stats as.formula cov formula model.matrix na.exclude na.omit predict resid simulate
bootstrap.lmerMod <- function(model, .f, type, B, resample, reb_type, linked){
  switch(type,
         parametric = parametric_bootstrap.lmerMod(model, .f, B, type = type),
         residual = resid_bootstrap.lmerMod(model, .f, B, type = type, linked = FALSE),
         case = case_bootstrap.lmerMod(model, .f, B, resample, type = type),
         cgr = cgr_bootstrap.lmerMod(model, .f, B, type = type),
         reb = reb_bootstrap.lmerMod(model, .f, B, reb_type = 0))
  # TODO: need to be able to save results
}


#' @rdname parametric_bootstrap
#' @export
parametric_bootstrap.lmerMod <- function(model, .f, B, type){
  .f <- match.fun(.f)
  
  model.fixef <- lme4::fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)
  
  # refit here
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  return(.bootstrap.completion(model, tstar, B, .f, type))
}


#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lmerMod <- function(model, .f, B, type, linked = FALSE){
  .f <- match.fun(.f)
  
  tstar <- purrr::map(1:B, function(x) .resample.resids(model, .f, linked = linked))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type)
  return(RES)
}


#' @rdname case_bootstrap
#' @export
case_bootstrap.lmerMod <- function(model, .f, B, resample, type){
  
  data <- model@frame
  # data$.id <- seq_len(nrow(data))
  clusters <- c(rev(names(lme4::getME(model, "flist"))), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  # rep.data <- purrr::map(integer(B), function(x) .cases.resamp(model = model, dat = data, cluster = clusters, resample = resample))
  tstar <- purrr::map(integer(B), function(x) .cases.resamp(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type)
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
.cases.resamp <- function(model, .f, dat, cluster, resample) {
  # exit early for trivial data
  if(nrow(dat) == 1 || all(resample==FALSE))
    return(dat)
  
  # ver <- as.numeric_version(packageVersion("dplyr"))
  res <- dat
  
  for(i in 1:length(cluster)) {
    if(i==1 & resample[i]) {
      dots <- as.name(cluster[1])
      grouped <- dplyr::group_by(res, !!dots)
      g_rows <- dplyr::group_rows(grouped)
      # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- as.name(cluster[-i])
        grouped <- dplyr::group_by(res, .dots = !!dots) 
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- as.name(cluster[i])
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- purrr::map_dfr(res, function(df) { # ldply to purrr map from list to df
            grouped <- dplyr::group_by(df, .dots = !!dots)
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
  
  if(class(model) == "lmerMod"){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    tstar <- .f(lme4::lmer(formula = form, data = res, REML = reml)) 
    return(tstar)
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(class(model) == "lme"){
    tstar <- tryCatch(.f(updated.model(model = model, new.data = res)),  
                      error = function(e) e)
    if(inherits(tstar, "error")) {
      structure(rep(NA, length(.f(model))), fail.msgs = tstar$message)
    } else{
      tstar
    }
    return(tstar)
  } else{
    stop("model class must be either 'lme' or 'lmerMod'")
  }
}


# .cases.resamp <- function (model, extra_step){
#   # Draw sample of size J from level-2 units
#   model.split <- split(x=model@frame, f=model@flist)
#   model.split.samp <- sample(x=model.split, size = length(model.split), resample = TRUE)
#   # For each sample, draw a sample of the cases from the level-2 unit
#   
#   if(extra_step == TRUE){
#     model.resamp <- lapply(model.split.samp,
#                            .f = function(x) {
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

# data parameter needs to change, will be result of .cases.resamp
# .cases.completion <- function(model, data = res, B, .f){
#   t0 <- .f(model)
#   
#   # moved to .cases.resamp
#   # # Refit the model and apply '.f' to it using map
#   
#   # form <- model@call$formula
#   # reml <- lme4::isREML(model)
#   # tstar <- purrr::map(data, function(x) {
#   #   .f(lme4::lmer(formula = form, data = x, REML = reml))
#   # })
#   
#   tstar <- res$tstar
#   tstar <- do.call("cbind", tstar) # Can these be nested?
#   rownames(tstar) <- names(t0)
#   
#   RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
#                         seed = .Random.seed, statistic = .f,
#                         sim = "case", call = match.call()),
#                    class = "boot")
#   
#   return(RES)
# }


#' @rdname cgr_bootstrap
#' @export
cgr_bootstrap.lmerMod <- function(model, .f, B, type){
  .f <- match.fun(.f)
  
  tstar <- as.data.frame(replicate(n = B, .resample.cgr(model = model, .f)))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type)
  return(RES)
}


#' @rdname reb_bootstrap
#' @export
reb_bootstrap.lmerMod <- function(model, .f, B, reb_type = 0){
  
  if(lme4::getME(object = model, name = "n_rfacs") > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  .f <- match.fun(.f)
  
  # changed this from ystar to tstar
  tstar <- as.data.frame(replicate(n = B, .resample.reb(model = model, .f = .f, reb_type = reb_type)))
  
  if(reb_type != 2) t0 <- .f(model)
  
  # Refit the model and apply '.f' to it using map
  #   tstar <- lapply(ystar[1,], function(x) {
  #     .f(refit(object = model, newresp = x))
  #   })
  
  if(reb_type == 2){
    fe.0 <- lme4::fixef(model)
    vc.0 <- bdiag(lme4::VarCorr(model))
    t0 <- c(beta = fe.0, sigma = c(diag(vc.0), lme4::getME(model, "sigma")^2))
    tstar <- purrr::map(tstar, function(x) { # changed ystar map to tstar map
      m <- lme4::refit(object = model, newresp = x)
      vc <- as.data.frame(lme4::VarCorr(m))
      list(fixef = lme4::fixef(m), varcomp = vc$vcov[is.na(vc$var2)])
    })
    
    vcs <- purrr::map(tstar, function(x) x$varcomp)
    Sb <- log( do.call("rbind", vcs) )
    #     fes <- lapply(tstar, function(x) x$fixef)
    
    Mb <- matrix(rep(colMeans(Sb), times = B), nrow = B, byrow = TRUE)
    CovSb <- cov(Sb)
    SdSb <- sqrt(diag(CovSb))
    
    Db <- matrix(rep(SdSb, times = B), nrow = B, byrow = TRUE)
    
    EW <- eigen(solve(CovSb), symmetric = T)
    Whalf <- EW$vectors %*% diag(sqrt(EW$values))
    
    Sbmod <- (Sb - Mb) %*% Whalf
    Sbmod <- Sbmod * Db # elementwise not a type 
    Lb <- exp(Mb + Sbmod)
    
    # tstar <- purrr::map(tstar, unlist)
  } else{
    Lb <- NULL
    # tstar <- purrr::map(ystar, function(x) {
    #   .f(lme4::refit(object = model, newresp = x))
    # })
  }
  
  # tstar <- do.call("cbind", tstar) # Can these be nested?
  #   rownames(tstar) <- names(.f(model))
  
  if(reb_type == 2) {
    idx <- 1:length(fe.0)
    fe.star <- tstar[idx,] 
    fe.adj <- sweep(fe.star, MARGIN = 1, STATS = fe.0 - rowMeans(fe.star, na.rm = TRUE), .f = "+")
    
    vc.star <- tstar[-idx,] 
    vc.adj <- sweep(vc.star, MARGIN = 1, STATS = t0[-idx] / rowMeans(vc.star, na.rm = TRUE), .f = "*")
    
    tstar <- rbind(fe.adj, vc.adj)
    
    .f <- function(.) {
      c(beta = lme4::fixef(.), sigma =c(diag(bdiag(lme4::VarCorr(.))), lme4::getME(., "sigma")^2))
    }
  }
  
  replicates <- as.data.frame(t(tstar))
  observed <- t0
  rep.mean <- colMeans(replicates)
  se <- unlist(purrr::map(replicates, sd))
  bias <- rep.mean - observed
  
  stats <- data.frame(observed, rep.mean, se, bias)
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = model@frame,
                        seed = .Random.seed, reb_type = reb_type, call = match.call()),
                   class = "lmeresamp")
  
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
#' This function combines \code{bstar} and \code{zstar} to create {Zbstar} using a map statement
#'
#' @param bstar A list of matrices bstar
#' @param zstar A list of matrices zstar
#'
#' @return matrix
#' @keywords internal
#' @noRd
.Zbstar.combine <- function(bstar, zstar){
  purrr::map(1:length(bstar), function(i){
    Matrix::t(zstar[[i]]) %*% bstar[[i]]
  })
}


#' @title Bootstrap Completion
#'
#' @description
#' Finishes the bootstrap process and makes the output readable.
#'
#' @details
#' This function is given \code{model, tstar, B, .f} and uses them to complete
#' the bootstrap process. They are then structured into a list for output and returned.
#'
#' @param tstar The tstar being passed in
#' @inheritParams bootstrap
#'
#' @return list
#' @keywords internal
#' @noRd
.bootstrap.completion <- function(model, tstar, B, .f, type = type){
  t0 <- .f(model)
  
  nsim <- length(tstar)
  tstar <- do.call("cbind", tstar) # Can these be nested?
  row.names(tstar) <- names(t0)
  
  if((nfail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", nfail, "/", nsim, ")")
    fail.msgs <- purrr::map_chr(tstar[bad.runs], .f = attr, FUN.VALUE = character(1),
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
                        stats = stats, R = B, data = model@frame,
                        seed = .Random.seed, type = type, call = match.call()), 
                   class = "lmeresamp")
  
  attr(RES,"bootFail") <- nfail
  attr(RES,"boot.fail.msgs") <- fail.msgs
  return(RES)
}

#' CGR resampling procedures
#' @keywords internal
#' @noRd
.resample.cgr <- function(model, .f){
  model.ranef <- lme4::ranef(model)
  
  # Extract residuals
  model.resid <- resid(model)
  
  # Higher levels
  Uhat.list <- purrr::map(seq_along(model.ranef),
                          .f = function(i) {
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
  ustar <- purrr::map(Uhat.list,
                      .f = function(df) {
                        index <- sample(x = seq_len(nrow(df)), size = nrow(df), replace = TRUE)
                        return(df[index,])
                      })
  
  # Structure u*
  if(level.num == 1){
    if(is.data.frame(ustar[[1]])){
      ustar <- purrr::map(ustar, .f = function(x) as.list(x))[[1]] 
    }
    names(ustar) <- names(Z)
  } else {
    ustar <- purrr::map(ustar, .f = function(x) as.data.frame(x))
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
  
  # Refit the model and apply '.f' to it using map
  # changed x to y.star for function and newresp
  tstar <- .f(lme4::refit(object = model, newresp = y.star))
}

#' Resampling residuals from mixed models
#' @param linkeded Specifies the use of marginal residual linkeding
#' @keywords internal
#' @noRd
.resample.resids <- function(model, .f, linked = linked) {
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  # centered
  model.ranef <- purrr::map(lme4::ranef(model), .f = scale, scale = FALSE)
  
  # Extract residuals
  # centered
  model.resid <- scale(resid(model), scale = FALSE)
  
  # Extract Z design matrix
  Z <- lme4::getME(object = model, name = "Ztlist")
  
  bstar <- purrr::map(model.ranef,
                      .f = function(x) {
                        J <- nrow(x)
                        
                        # Sample of b*
                        bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                        bstar <- x[bstar.index,]
                        return(bstar)
                      })
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  if(level.num == 1){
    if(!is.numeric(bstar[[1]])) bstar <- purrr::map(bstar, .f = function(x) as.list(x))[[1]]
    names(bstar) <- names(Z)
  } else {
    bstar <- purrr::map_dfr(bstar, .f = function(x) as.data.frame(x))
    bstar <- do.call(c, bstar)
    names(bstar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  if(linked == FALSE){
    # Resample residuals
    estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
    
    # Combine function
    y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
    
    # Refit the model and apply '.f' to it using map
    
    ## .f() no mapping
    tstar <- .f(lme4::refit(object = model, newresp = y.star))
  } else{
    # linked
    model.mresid <- lme4::getME(model, "y") - predict(model, re.form = ~0)
    model.mresid.cent <- scale(model.mresid, scale = FALSE)
    
    # Resample residuals
    mresid.star <- sample(x = model.mresid.cent, size = length(model.mresid.cent), replace = TRUE)
    
    # Combine function
    y.star <- as.numeric(Xbeta + mresid.star)
    
    # Refit the model and apply '.f' to it using map
    
    ## .f() no mapping
    tstar <- .f(lme4::refit(object = model, newresp = y.star))
  }
}

#' REB resampling procedures 
#' 
#' @param reb_type Specifies the inclusion of REB/1
#' @inheritParams bootstrap
#' @import Matrix
#' @keywords internal
#' @noRd
.resample.reb <- function(model, .f, reb_type){
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
  levs <- purrr::map(fl <- model@flist, levels)
  asgn <- attr(fl, "assign")
  cnms <- model@cnms
  nc <- vapply(cnms, length, 1L)
  nb <- nc * (nl <- vapply(levs, length, 1L)[asgn])
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(u, nbseq)
  for (i in seq_along(u)){
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  }
  names(u) <- names(cnms)
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  if(reb_type == 1){
    if(level.num > 1) stop("reb_type = 1 is not yet implemented for higher order models")
    # Calculations
    Uhat <- purrr::map(u, function(x){
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
  ustar <- purrr::map(Uhat,
                      .f = function(x) {
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
    ustar <- purrr::map(ustar, .f = function(x) as.list(x))[[1]]
  } else {
    ustar <- purrr::map_dfr(ustar, .f = function(x) as.data.frame(x))
    ustar <- do.call(c, ustar)
  }
  
  names(ustar) <- names(Ztlist) 
  ustar.df <- as.data.frame(ustar)
  
  #   if(level.num == 1){
  #     Uhat.list <- lapply(Uhat.list, .f = function(x) as.list(x))[[1]]
  #     names(Uhat.list) <- names(Ztlist)
  #   } else {
  #     Uhat.list <- sapply(Uhat.list, .f = function(x) as.list(x))
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
  
  # Refit the model and apply '.f' to it using map
  tstar <- .f(lme4::refit(object = model, newresp = y.star))
}

