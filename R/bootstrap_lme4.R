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



