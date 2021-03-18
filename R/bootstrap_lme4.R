#' @rdname bootstrap
#' @export
#' @importFrom stats as.formula cov formula model.matrix na.exclude 
#' na.omit predict resid simulate sd confint quantile
bootstrap.merMod <- function(model, .f, type, B, resample, reb_type, hccme, aux.dist){
  switch(type,
         parametric = parametric_bootstrap.merMod(model, .f, B),
         residual = resid_bootstrap.merMod(model, .f, B),
         case = case_bootstrap.merMod(model, .f, B, resample),
         reb = reb_bootstrap.lmerMod(model, .f, B, reb_type),
         wild = wild_bootstrap.lmerMod(model, .f, B, hccme, aux.dist))
}


#' @rdname parametric_bootstrap
#' @export
parametric_bootstrap.merMod <- function(model, .f, B){
  .f <- match.fun(.f)
  
  # model.fixef <- lme4::fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)
  
  # refit here
  tstar <- purrr::map(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  return(.bootstrap.completion(model, tstar, B, .f, type = "parametric"))
}


#' @rdname case_bootstrap
#' @export
case_bootstrap.merMod <- function(model, .f, B, resample){
  
  data <- model@frame
  flist <- lme4::getME(model, "flist")
  re_names <- names(flist)
  clusters <- c(rev(re_names), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. 
         Please specify whether to resample the data at each level of grouping,
         including at the observation level.")
  
  if(!all(re_names %in% colnames(data))) {
    missing_re <- setdiff(re_names, colnames(data))
    data <- dplyr::bind_cols(data, flist[missing_re])
  }
  
  # rep.data <- purrr::map(integer(B), function(x) .cases.resamp(model = model, dat = data, cluster = clusters, resample = resample))
  tstar <- purrr::map(integer(B), function(x) .resample_refit.cases(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type = "case")
  return(RES)
}



#' @rdname resid_bootstrap
#' @export
resid_bootstrap.merMod <- function(model, .f, B){
  
  .f <- match.fun(.f)
  glmm <- lme4::isGLMM(model)
  
  setup <- .setup(model, type = "residual")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.cgr(
        glmm = glmm,
        b = setup$b, 
        e = setup$e, 
        level.num = setup$level.num, 
        Ztlist = setup$Ztlist, 
        Xbeta = setup$Xbeta,
        vclist = setup$vclist,
        sig0 = setup$sig0,
        invlink = ifelse(glmm, model@resp$family$linkinv, NULL)
      )
    )
  )
  
  if(glmm){
    fam <- stats::family(model)
    wts <- stats::weights(model)
    
    # simulate y
    simfun <- simfunList[[fam$family]]
    ystar <- purrr::map(ystar,
      ~simfun(model, nsim = 1, ftd = .x, wts = wts)
    )
      
  }
  
  tstar <- purrr::map(ystar, function(x) {
    .f(lme4::refit(object = model, newresp = x))
  })
  
  .bootstrap.completion(model, tstar, B, .f, type = "residual")
}


#' @rdname wild_bootstrap
#' @export
wild_bootstrap.lmerMod <- function(model, .f, B, hccme = c("hc2", "hc3"), 
                                   aux.dist = c("f1", "f2")){
  
  .f <- match.fun(.f)
  hccme <- match.arg(hccme)
  aux.dist <- match.arg(aux.dist)
  
  setup <- .setup(model, type = "wild")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.wild(
        Xbeta = setup$Xbeta, 
        mresid = setup$mresid, 
        .hatvalues = setup$.hatvalues, 
        hccme = hccme, 
        aux.dist = aux.dist,
        n.lev = setup$n.lev,
        flist = setup$flist
      )
    )
  )
  
  tstar <- purrr::map(ystar, function(y) {
    .f(lme4::refit(object = model, newresp = y))
  })
  
  .bootstrap.completion(model, tstar, B, .f, type = "wild")
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
  if(reb_type == 2) .f <- extract_parameters.merMod
  
  tstar <- purrr::map(y.star, function(x) {
    .f(lme4::refit(object = model, newresp = x))
    })
  
  # Extract original statistics
  t0 <- .f(model)
  
  # Postprocessing for REB/2
  if(reb_type == 2) 
    tstar <- .postprocess.reb2(t0, tstar, nbeta = length(lme4::getME(model, "beta")), B = B)
  
  # Format for return
  .bootstrap.completion(model, tstar, B, .f, type = paste("reb", reb_type, sep = ""))
}



