#' Case resampler for mixed models
#' @keywords internal
#' @noRd
.resamp.cases <- function(dat, cluster, resample) {
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
  
  res
}


.resample_refit.cases <- function(model, .f, dat, cluster, resample){
  resamp_data <- .resamp.cases(dat, cluster, resample)
  
  if(class(model) == "lmerMod"){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    tstar <- .f(lme4::lmer(formula = form, data = resamp_data, REML = reml)) 
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(class(model) == "lme"){
    tstar <- tryCatch(.f(updated.model(model = model, new.data = resamp_data)),  
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



#' CGR resampling from lmerMod objects
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

#' Resampling residuals from lmerMod objects
#' @keywords internal
#' @noRd
.resample.resids <- function(b, e, level.num, Ztlist, Xbeta) {
  bstar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  
  if(level.num == 1){
    if(!is.numeric(bstar[[1]])) bstar <- purrr::map(bstar, .f = as.list)[[1]]
    names(bstar) <- names(Ztlist)
  } else {
    # bstar <- purrr::map_dfr(bstar, .f = as.data.frame)
    # bstar <- do.call(c, bstar)
    # names(bstar) <- names(Ztlist)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # Resample residuals
  estar <- sample(x = e, size = length(e), replace = TRUE)
  
  # Calculate boostrap responses
  as.numeric(Xbeta + Zbstar.sum + estar)
}


#' REB resampling procedure for lmerMod objects
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


#' Resampling residuals from lme objects
#' @keywords internal
#' @noRd
.resample.resids.lme <- function(b, e, Xbeta, Zlist){
  bstar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  estar <- sample(e, size = length(e), replace = TRUE)
  
  Zbstar.sum <- .Zbstar.combine.lme(bstar, Zlist)
  
  # Combine function
  as.numeric(Xbeta + Zbstar.sum + estar)
}

#' CGR resampler for lme objects
#' @keywords internal
#' @noRd
.resample.cgr.lme <- function(b, e, Xbeta, Zlist, vclist, sig0){
  # Resample resids
  estar <- sample(x = e, size = length(e), replace = TRUE)
  ehat <- scale_center_e(estar, sigma = sig0)
  
  # Resample Uhat
  ustar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  ustar <- purrr::map2(ustar, vclist, scale_center_ranef)
  
  Zbstar.sum <- .Zbstar.combine.lme(ustar, Zlist)
  
  # Combine function
  as.numeric(Xbeta + Zbstar.sum + estar)
}


#' REB resampler for lme objects
#' @keywords internal
#' @noRd
.resample.reb.lme <- function(Xbeta, Zlist, Uhat, estar.vec, flist, levs){
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
  # ustar <- ustar[[1]]
  # 
  # names(ustar) <- names(Zlist) 
  # ustar.df <- as.data.frame(ustar)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine.lme(bstar = ustar, Zlist = Zlist)
  # Zbstar.sum <- Reduce("+", Zbstar)
  
  # ystar
  as.numeric(Xbeta + Zbstar + estar)
}