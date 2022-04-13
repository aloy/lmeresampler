# bootstrap_resample.merMod <- function(model, type, B, resample, reb_type, hccme, 
#                                aux.dist, orig_data = NULL) {
#   
#   if(type == "parametric") {
#     ystar <- simulate(model, nsim = B, na.action = na.exclude)
#   }
#   
#   if(type == "cases") {
#     prep <- prep_cases.merMod(model = model, resample = resample, orig_data = orig_data)
#   }
#   
#   if(type == "residual") {
#     
#   }
#   
#   
# }



#' Case resampler for mixed models
#' @keywords internal
#' @noRd
.resamp.cases <- function(dat, cluster, resample) {
  res <- dat
  
  for(i in seq_along(cluster)) {
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


.resample_refit.cases <- function(model, .f, dat, cluster, resample, .refit){
  resamp_data <- .resamp.cases(dat, cluster, resample)
  error <- NULL
  
  if(!.refit) return(resamp_data)
  
  if(inherits(model, "lmerMod")){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    f1 <- factory(
      function(form, resamp_data, reml) 
        .f(lme4::lmer(formula = form, data = resamp_data, REML = reml))
      )
    tstar <- f1(form, resamp_data, reml)
    
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(inherits(model, "lme")){
    tstar <- updated.model(model = model, new.data = resamp_data)  
    tstar <- .f(tstar)
  } else if(inherits(model, "glmerMod")) {
    form <- update(model@call$formula, y ~ .)
    y_idx <- colnames(resamp_data) == getResponseFromFormula(model)
    if(sum(y_idx) > 0) {
      colnames(resamp_data)[y_idx] <- "y"
    } else {
      colnames(resamp_data)[1] <- "y"
    }
    fam  <- family(model)
    
    f1 <- factory(
      function(form, resamp_data, fam) 
        .f(lme4::glmer(formula = form, data = resamp_data, family = fam))
    )
    tstar <- f1(form, resamp_data, fam)
    
  } else{
    stop("model class must be one of 'lme', 'lmerMod', or 'glmerMod'")
  }
  tstar
}



#' CGR resampling from lmerMod objects
#' @keywords internal
#' @noRd
.resample.cgr <- function(glmm, b, e, level.num, Ztlist, Xbeta, vclist, sig0, invlink){
  
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
  
  if(glmm) {
    eta <- as.numeric(Xbeta + Zbstar.sum)
    ystar  <- invlink(eta) # not really ystar, need inv. link
  } else{
  # Resample resids
  estar <- sample(x = e, size = length(e), replace = TRUE)
  ehat <- scale_center_e(estar, sigma = sig0)
  
  # Calc. bootstrap y
  ystar <- as.numeric(Xbeta + Zbstar.sum + estar)
  }
  
  ystar
}



#' Wild bootstrap resampling from lmerMod and lme objects
#' @keywords internal
#' @noRd
.resample.wild <- function(Xbeta, mresid, .hatvalues, hccme, aux.dist, 
                           flist, n.lev){
  
  # Sample from auxillary distribution
  # Mammen distribution
  if(aux.dist == "mammen") { 
    prob <- (sqrt(5) + 1) / (2 * sqrt(5))
    w <- sample(
      c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2), 
      size = n.lev,
      prob = c(prob, 1 - prob),
      replace = TRUE
    )
  } 
  
  # Rademacher distribution
  if(aux.dist == "rademacher") {
    w <- sample(c(1, -1), size = n.lev, replace = TRUE)
  }
  
  # Standard normal
  if(aux.dist == "norm") {
    w <- rnorm(n = n.lev)
  }
  
  # Webb's 6-point distribution
  if(aux.dist == "webb") {
    w <- sample(
      c(sqrt(3/2), 1, sqrt(1/2), -sqrt(3/2), -1, -sqrt(1/2)),
      size = n.lev,
      replace = TRUE
    )
  }
  
  # Recentered Gamma(4, scale = 1/2) -- Liu (1988)
  if(aux.dist == "gamma") {
    w <- rgamma(n = n.lev, shape = 4, scale = 1/2) - 2
  }
  
  # Calc. bootstrap y
  if(hccme == "hc2") v <- (1 / sqrt(1 - .hatvalues)) * mresid
  if(hccme == "hc3") v <- (1 / (1 - .hatvalues)) * mresid
    
    as.numeric(Xbeta + v * w[flist])
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
