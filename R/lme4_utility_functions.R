scale_center_ranef <- function(b, vc){
  u <- scale(b, scale = FALSE)
  S <- (t(u) %*% u) / nrow(u)
  R <- vc
  
  Ls <- t(chol(S, pivot = TRUE))
  Lr <- t(chol(R, pivot = TRUE))
  A <- t(Lr %*% solve(Ls))
  
  Uhat <- as.matrix(u %*% A)
  data.frame(Uhat)
}

scale_center_e <- function(x, sigma) {
  estar <- sigma * x %*% ((t(x) %*% x) / length(x))^(-1/2)
  scale(estar, scale = FALSE)
}


#' @importFrom stats sigma
#' @export
extract_parameters.merMod <- function(model) {
  sig.e <- stats::sigma(model)
  vc <- as.data.frame(lme4::VarCorr(model))
  
  c(
    beta = lme4::getME(model, "beta"), 
    vc = vc$vcov[is.na(vc$var2)]
  )
}


#' Organize ranef vector into lists
#' 
#' @param b ranef estimates
#' @param fl flist from lmerMod object, a list of the grouping 
#'   variables (factors) involved in the random effect terms
#' @param levs unique levels of each factor in fl
#' @param cnms a list component names for each ranef
#' @keywords internal
#' @noRd
arrange_ranefs.lmerMod <- function(b, fl, levs, cnms){
  asgn <- attr(fl, "assign")
  nc <- vapply(cnms, length, 1L)
  nb <- nc * (nl <- vapply(levs, length, 1L)[asgn])
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(b, nbseq)
  for (i in seq_along(u)){
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  }
  names(u) <- names(cnms)
  u
}


arrange_ranefs.lme <- function(b, fl, levs, cnms){
  asgn <- seq_along(fl)
  nc <- purrr::map_int(cnms, length) #map_int()
  nb <- nc * (nl <- purrr::map_int(levs, length))
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(b, nbseq)
  for (i in seq_along(u)){
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  }
  names(u) <- names(cnms)
  u
}


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
    Matrix::t(zstar[[i]]) %*% as.matrix(bstar[[i]])
  })
}



#' REB/2 post processing
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
  
  tstar <- dplyr::bind_cols(tstar)
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
  
  res <- as.data.frame(t(cbind(fe.adj, vc.adj)))
  res <- lapply(seq(nrow(res)), function(i) as.numeric(res[i, ]))
  names(res) <- names(t0)
  
  res
}

#' Refitting merMod with error catching
#' @param ystar bootstrapped responses
#' @param model fitted merMod object
#' @param .f function to calc bootstrap stats
#' @keywords internal
#' @noRd
#' @importFrom stats napredict
refit_merMod <- function(ystar, model, .f) {
  error <- NULL
  
  # Adjustment to respect na.action
  .na.act <- attr(model@frame, "na.action")
  ystar2 <- purrr::map(ystar, function(.y) {
    attr(.y, "na.action") <- .na.act
    .y
  })
  
  
  f1 <- factory(
    function(model, y) 
      .f(lme4::refit(object = model, newresp = y))
  )
  stats <- purrr::map(ystar2, function(.y) f1(model, .y))

  
  list(tstar = stats, warnings = collect_warnings(stats))
}

#' @importFrom stats terms getCall
getResponseFromFormula <- function(model) {
  if (attr(terms(model), which = 'response'))
    all.vars(getCall(model)$formula)[1]
  else
    NULL
}

refit_to_newdf <- function(model, newdata, .f) {
  if(inherits(model, "lmerMod")){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    f1 <- factory(
      function(form, newdata, reml) 
        .f(lme4::lmer(formula = form, data = newdata, REML = reml))
    )
    tstar <- f1(form, newdata, reml)
    
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(inherits(model, "lme")){
    tstar <- updated.model(model = model, new.data = newdata)  
    tstar <- .f(tstar)
  } else if(inherits(model, "glmerMod")) {
    form <- update(model@call$formula, y ~ .)
    y_idx <- colnames(newdata) == getResponseFromFormula(model)
    if(sum(y_idx) > 0) {
      colnames(newdata)[y_idx] <- "y"
    } else {
      colnames(newdata)[1] <- "y"
    }
    fam  <- family(model)
    
    f1 <- factory(
      function(form, newdata, fam) 
        .f(lme4::glmer(formula = form, data = newdata, family = fam))
    )
    tstar <- f1(form, newdata, fam)
    
  } else{
    stop("model class must be one of 'lme', 'lmerMod', or 'glmerMod'")
  }
  tstar
}


#' Set up for case resampling in lme4
# prep_cases.merMod <- function(model, resample, orig_data) {
#   if(!is.null(orig_data)){
#     data <- orig_data
#   }else{
#     data <- model@frame
#   }
#   
#   flist <- lme4::getME(model, "flist")
#   re_names <- names(flist)
#   clusters <- c(rev(re_names), ".id")
#   
#   if(length(clusters) != length(resample))
#     stop("'resample' is not the same length as the number of grouping variables. 
#          Please specify whether to resample the data at each level of grouping,
#          including at the observation level.")
#   
#   if(!all(re_names %in% colnames(data))) {
#     missing_re <- setdiff(re_names, colnames(data))
#     data <- dplyr::bind_cols(data, flist[missing_re])
#   }
#   
#   return(list(data, clusters))
# }
