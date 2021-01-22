#' @rdname bootstrap
#' @export
#' @importFrom stats as.formula cov formula model.matrix na.exclude 
#' na.omit predict resid simulate sd confint quantile
bootstrap.lmerMod <- function(model, .f, type, B, resample, reb_type){
  switch(type,
         parametric = parametric_bootstrap.lmerMod(model, .f, B),
         residual = resid_bootstrap.lmerMod(model, .f, B),
         case = case_bootstrap.lmerMod(model, .f, B, resample),
         cgr = cgr_bootstrap.lmerMod(model, .f, B),
         reb = reb_bootstrap.lmerMod(model, .f, B, reb_type))
  # TODO: need to be able to save results
}


#' @rdname parametric_bootstrap
#' @export
parametric_bootstrap.lmerMod <- function(model, .f, B){
  .f <- match.fun(.f)
  
  # model.fixef <- lme4::fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)
  
  # refit here
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  return(.bootstrap.completion(model, tstar, B, .f, type = "parametric"))
}


#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lmerMod <- function(model, .f, B){
  
  .f <- match.fun(.f)
  
  setup <- .setup(model, type = "resids")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.resids(
        b = setup$b, 
        e = setup$e, 
        level.num = setup$level.num, 
        Ztlist = setup$Ztlist, 
        Xbeta = setup$Xbeta
      )
    )
  )
  
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  
  .bootstrap.completion(model, tstar, B, .f, type = "residual")
}


#' @rdname case_bootstrap
#' @export
case_bootstrap.lmerMod <- function(model, .f, B, resample){
  
  data <- model@frame
  # data$.id <- seq_len(nrow(data))
  clusters <- c(rev(names(lme4::getME(model, "flist"))), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  # rep.data <- purrr::map(integer(B), function(x) .cases.resamp(model = model, dat = data, cluster = clusters, resample = resample))
  tstar <- purrr::map(integer(B), function(x) .resamp.cases(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type = "case")
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
.resamp.cases <- function(model, .f, dat, cluster, resample) {
  
  res <- dat
  
  for(i in 1:length(cluster)) {
    if(i==1 & resample[i]) {
      dots <- as.name(cluster[1])
      grouped <- dplyr::group_by(res, !!dots)
      g_rows <- dplyr::group_rows(grouped)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- as.name(cluster[-i])
        grouped <- dplyr::group_by(res, !!dots) 
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- as.name(cluster[i])
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- purrr::map_dfr(res, function(df) { # ldply to purrr map from list to df
            grouped <- dplyr::group_by(df, !!dots)
            g_rows <- dplyr::group_rows(grouped)
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
cgr_bootstrap.lmerMod <- function(model, .f, B){
  
  .f <- match.fun(.f)
  
  setup <- .setup(model, type = "cgr")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.cgr(
        b = setup$b, 
        e = setup$e, 
        level.num = setup$level.num, 
        Ztlist = setup$Ztlist, 
        Xbeta = setup$Xbeta,
        vclist = setup$vclist,
        sig0 = setup$sig0
      )
    )
  )
  
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  
  .bootstrap.completion(model, tstar, B, .f, type = "cgr")
}


#' @rdname reb_bootstrap
#' @export
reb_bootstrap.lmerMod <- function(model, .f, B, reb_type){
  
  if(missing(reb_type)){
    reb_type <- 0
    warning("'reb_type' unspecified, performing REB 0 bootstrap")
  }
  
  if(lme4::getME(object = model, name = "n_rfacs") > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  .f <- match.fun(.f)
  
  # Set up for bootstrapping
  setup <- .setup(model, type = "reb", reb_type = reb_type)
  
  # Generate bootstrap responses
  y.star <- replicate(
    n = B, 
    .resample.reb(
      Xbeta = setup$Xbeta, 
      Ztlist = setup$Ztlist, 
      Uhat = setup$b, 
      estar.vec = as.numeric(setup$e), 
      flist = setup$flist, 
      levs = setup$levs
    )
  )
  
  y.star <- as.data.frame(y.star)
  
  # Extract bootstrap statistics
  if(reb_type == 2) .f <- extract_parameters.lmerMod
  
  tstar <- purrr::map_dfc(y.star, function(x) {
    .f(lme4::refit(object = model, newresp = x))
    })
  
  # Extract original statistics
  t0 <- .f(model)
  
  # Postprocessing for REB/2
  if(reb_type == 2) 
    tstar <- .postprocess.reb2(t0, tstar, nbeta = length(getME(model, "beta")), B = B)
  
  # Format for return
  .bootstrap.completion(model, tstar, B, .f, type = paste("reb", reb_type, sep = ""))
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
  
  if(class(model) == "lmerMod") {
    data = model@frame
  } else if (class(model) == "lme") {
    data = model$data
  }
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, R = B, data = data,
                        seed = .Random.seed, type = type, call = match.call()), 
                   class = "lmeresamp")
  
  attr(RES,"bootFail") <- nfail
  attr(RES,"boot.fail.msgs") <- fail.msgs
  return(RES)
}



#' CGR resampling procedures
#' @keywords internal
#' @noRd
.resample.cgr <- function(b, e, level.num, Ztlist, Xbeta, vclist, sig0){
  
  # Resample resids
  estar <- sample(x = e, size = length(e), replace = TRUE)
  ehat <- scale_center_e(estar, sigma = sig0)
  
  # Resample Uhat
  ustar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  ustar <- purrr::map2(ustar, vclist, scale_center_ranef)
  
  # Structure u*
  if(level.num == 1){
    if(is.data.frame(ustar[[1]])){
      ustar <- purrr::map(ustar, .f = as.list)[[1]]
    }
    names(ustar) <- names(Ztlist)
  } else {
    ustar <- purrr::map(ustar, .f = as.data.frame)
    ustar <- do.call(c, ustar)
    names(ustar) <- names(Ztlist)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Calc. bootstrap y
  as.numeric(Xbeta + Zbstar.sum + estar)
}

#' Resampling residuals from mixed models
#' @keywords internal
#' @noRd
.resample.resids <- function(b, e, level.num, Ztlist, Xbeta) {
  bstar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  
  if(level.num == 1){
    if(!is.numeric(bstar[[1]])) bstar <- purrr::map(bstar, .f = as.list)[[1]]
    names(bstar) <- names(Ztlist)
  } else {
    bstar <- purrr::map_dfr(bstar, .f = as.data.frame)
    bstar <- do.call(c, bstar)
    names(bstar) <- names(Ztlist)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = e, size = length(e), replace = TRUE)
  
  # Calculate boostrap responses
  as.numeric(Xbeta + Zbstar.sum + estar)
}




#' REB resampling procedure 
#' 
#' @param Xbeta marginal fitted values
#' @param Ztlist design matrix separated by variance
#' @param Uhat ranefs organized as Ztlist
#' @param estar.vec vector of level-1 residuals
#' @param flist a list of the grouping variables (factors) involved in the random effect terms
#' @param levs a list of levels of the grouping variables in flist
#' @inheritParams bootstrap
#' @import Matrix
#' @keywords internal
#' @noRd
.resample.reb <- function(Xbeta, Ztlist, Uhat, estar.vec, flist, levs){
  # For now, assume only a two-level model
  grps <- levs[[1]]
  units <- flist[[1]]
  resamp_u_ids <- sample(seq_along(grps), size = length(grps), replace = TRUE)
  resamp_e_grps <- sample(grps, size = length(grps), replace = TRUE)
  
  # resample uhats
  ustar <- purrr::map(Uhat, ~data.frame(.x[resamp_u_ids, ]))
  
  # Resample residuals, e
  estar <- numeric(length = length(units))
  for(i in seq_along(resamp_e_grps)) {
    target_units <- which(units == grps[i])
    donor_units <- which(units == resamp_e_grps[i])
    estar[target_units] <- sample(estar.vec, size = length(target_units), replace = TRUE)
  }

  # since only working with 2-levels models now
  ustar <- ustar[[1]]
  
  names(ustar) <- names(Ztlist) 
  ustar.df <- as.data.frame(ustar)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar.df, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # ystar
  as.numeric(Xbeta + Zbstar.sum + estar)
}


#' REB resampling procedure 
#' 
#' @param t0 stats extracted from original model
#' @param tstar bootstrap statistics
#' @param nbeta number of fixed effects parameters, \code{length(getME(model, "beta"))}
#' @keywords internal
#' @noRd
.postprocess.reb2 <- function(t0, tstar, nbeta, B){
  nparams <- length(t0)
  fe0 <- t0[1:nbeta]
  vc0 <- t0[(nbeta+1):nparams]
  
  fe.star <- t(tstar[1:nbeta,])
  vc.star <- t(tstar[(nbeta+1):nparams,])
  
  # Rescaling
  Sb <- log(vc.star)
  CovSb <- cov(Sb)
  SdSb <- sqrt(diag(CovSb))
  
  EW <- eigen(solve(CovSb), symmetric = T)
  Whalf <- EW$vectors %*% diag(sqrt(EW$values))
  
  Db <- matrix(rep(SdSb, times = B), nrow = B, byrow = TRUE)
  Mb <- matrix(rep(colMeans(Sb), times = B), nrow = B, byrow = TRUE)
  
  Sbmod <- (Sb - Mb) %*% Whalf
  Sbmod <- Sbmod * Db # elementwise not a type 
  Lb <- exp(Mb + Sbmod)
  
  # Bias corrections
  fe.adj <- sweep(fe.star, MARGIN = 2, STATS = fe0 - colMeans(fe.star, na.rm = TRUE), FUN = "+")
  vc.adj <- sweep(vc.star, MARGIN = 2, STATS = vc0 / colMeans(vc.star, na.rm = TRUE), FUN = "*")
  
  as.data.frame(t(cbind(fe.adj, vc.adj)))
}

