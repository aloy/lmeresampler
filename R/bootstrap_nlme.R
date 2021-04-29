#' @rdname bootstrap
#' @export
bootstrap.lme <- function(model, .f, type, B, resample, reb_type, hccme, aux.dist){
  switch(type,
         parametric = parametric_bootstrap.lme(model, .f, B),
         residual = resid_bootstrap.lme(model, .f, B),
         case = case_bootstrap.lme(model, .f, B, resample),
         reb = reb_bootstrap.lme(model, .f, B, reb_type),
         wild = wild_bootstrap.lme(model, .f, B, hccme, aux.dist))
}


#' @rdname parametric_bootstrap
#' @export
#' @importFrom nlmeU simulateY
parametric_bootstrap.lme <- function(model, .f, B){
  # Match function
  .f <- match.fun(.f)
  t0 <- .f(model)
  
  # Extract fixed effects
  # model.fixef <- nlme::fixef(model)
  
  ystar <- nlmeU::simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  ystar <- data.frame(ystar)
  
  refits <- refit_lme(ystar = ystar, model = model, .f = .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "parametric", warnings = refits$warnings)
}



#' @rdname case_bootstrap
#' @export
case_bootstrap.lme <- function(model, .f, B, resample){
  
  data <- model$data
  # data$.id <- seq_len(nrow(data))
  clusters <- c(names(model$groups), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. 
         Please specify whether to resample the data at each level of grouping.")
  
  refits <- purrr::map(integer(B), function(x) .resample_refit.cases(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  tstar <- purrr::map(refits, ~.x$value)
  warnings <- list(
    warning = lapply(refits, function(.x) unlist(.x$warning)$message),
    message = lapply(refits, function(.x) unlist(.x$message)$message),
    error = lapply(refits, function(.x) unlist(.x$error)$message)
  )
  
  .bootstrap.completion(model, tstar = tstar, B, .f, type = "case", warnings = warnings)
}



#' @rdname reb_bootstrap
#' @inheritParams bootstrap
#' @export
reb_bootstrap.lme <- function(model, .f, B, reb_type){
  
  if(missing(reb_type)){
    reb_type <- 0
    warning("'reb_type' unspecified, performing REB 0 bootstrap")
  }
  
  if(ncol(model$groups) > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  .f <- match.fun(.f)
  
  # Set up for bootstrapping
  setup <- .setup.lme(model, type = "reb", reb_type = reb_type)
  
  # Generate bootstrap responses
  y.star <- replicate(
    n = B, 
    .resample.reb.lme(
      Xbeta = setup$Xbeta, 
      Zlist = setup$Zlist, 
      Uhat = setup$b, 
      estar.vec = as.numeric(setup$e), 
      flist = setup$flist, 
      levs = setup$levs
    )
  )
  
  ystar <- as.data.frame(y.star)
  
  # Extract bootstrap statistics
  if(reb_type == 2) .f <- extract_parameters.lme
  
  refits <- refit_lme(ystar = ystar, model = model, .f = .f)
  tstar <- refits$tstar
  
  # Extract original statistics
  t0 <- .f(model)
  
  # Postprocessing for REB/2
  if(reb_type == 2) 
    tstar <- .postprocess.reb2(t0, tstar, nbeta = length(nlme::fixef(model)), B = B)
  
  # Format for return
  .bootstrap.completion(model, tstar = tstar, B, .f, 
                        type = paste("reb", reb_type, sep = ""), 
                        warnings = refits$warnings)
}





#' @rdname resid_bootstrap
#' @inheritParams bootstrap
#' @export
resid_bootstrap.lme <- function(model, .f, B){
  
  .f <- match.fun(.f)
  
  setup <- .setup.lme(model, type = "residual")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.cgr.lme(
        b = setup$b, 
        e = setup$e, 
        Xbeta = setup$Xbeta,
        Zlist = setup$Zlist,
        vclist = setup$vclist, 
        sig0 = setup$sig0
      )
    )
  )
  
  
  refits <- refit_lme(ystar = ystar, model = model, .f = .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "residual", warnings = refits$warnings)
}


#' @rdname wild_bootstrap
#' @inheritParams bootstrap
#' @export
wild_bootstrap.lme <- function(model, .f, B, hccme = c("hc2", "hc3"), 
                               aux.dist = c("f1", "f2")){
  
  .f <- match.fun(.f)
  hccme <- match.arg(hccme)
  aux.dist <- match.arg(aux.dist)
  
  setup <- .setup.lme(model, type = "wild")
  
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
  
  
  refits <- refit_lme(ystar = ystar, model = model, .f = .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "wild", warnings = refits$warnings)
}

