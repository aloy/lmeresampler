#' Setup for resampling
#' 
#' @inheritParams bootstrap
#' @keywords internal
#' @noRd
.setup <- function(model, type, reb_type = NULL){
  # Extract marginal means
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract ranef design matrix list
  Ztlist <- lme4::getME(object = model, name = "Ztlist")
  
  level.num <- lme4::getME(object = model, name = "n_rfacs")
  
  if(type == "cgr" || type == "reb" && reb_type == 1){
    sig0 <- sigma(model)
    vclist <- purrr::map(
      seq_along(b), 
      .f = ~bdiag(lme4::VarCorr(model)[[names(b)[.x]]])
    )
    names(vclist) <- names(b)
  }
  
  if(type == "reb") {
    Z <- lme4::getME(object = model, name = "Z")
    mresid <- lme4::getME(model, "y") - Xbeta
    
    # level-2 resid
    b <- solve(t(Z) %*% Z) %*% t(Z) %*% mresid # a single vector
    
    # level 1 resid
    e <- mresid - Z %*% b
    
    levs <- purrr::map(fl <- model@flist, levels)
    cnms <- getME(model, "cnms")
    b <- arrange_ranefs.lmerMod(b, fl, levs, cnms)
    
    if(reb_type == 1){
      if(level.num > 1) stop("reb_type = 1 is not yet implemented for higher order models")
      
      # Rescale u the residuals *prior* to resampling
      # so empirical variance is equal to estimated variance
      b <- purrr::map2(b, vclist, scale_center_ranef)
      e <- scale_center_e(e, sig0)
    } 
  } else{
    # Extract and center random effects
    b <- purrr::map(lme4::ranef(model), .f = scale, scale = FALSE)
    b <- purrr::map(b, as.data.frame)
    
    # Extract and center error terms
    e <- scale(resid(model), scale = FALSE)
    
    if(type == "cgr"){
      sig0 <- sigma(model)
      vclist <- purrr::map(
        seq_along(b), 
        .f = ~bdiag(lme4::VarCorr(model)[[names(b)[.x]]])
      )
      names(vclist) <- names(b)
    }
  }
  
  
  RES <- list(Xbeta = Xbeta, b = b, e = e, Ztlist = Ztlist)
  if(type == "reb") {
    RES <- append(RES, list(flist = fl, levs = levs))
  } else{
    RES <- append(RES, list(level.num = level.num))
    if(type == "cgr") RES <- append(RES, list(sig0 = sig0, vclist = vclist))
  }
  
  RES
}
