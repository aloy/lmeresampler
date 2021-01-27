#' @rdname bootstrap
#' @export
bootstrap.lme <- function(model, .f, type, B, resample, reb_type){
  switch(type,
         parametric = parametric_bootstrap.lme(model, .f, B),
         residual = resid_bootstrap.lme(model, .f, B),
         case = case_bootstrap.lme(model, .f, B, resample),
         cgr = cgr_bootstrap.lme(model, .f, B),
         reb = reb_bootstrap.lme(model, .f, B, reb_type))
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
  
  tstar <- purrr::map(ystar, function(y) {
    fit <- tryCatch(.f(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if(inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  })
  
  return(.bootstrap.completion(model, tstar, B, .f, type = "parametric"))
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
  
  tstar <- purrr::map(integer(B), function(x) .resamp.cases(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  return(.bootstrap.completion(model, tstar, B, .f, type = "case"))
}

#' @rdname resid_bootstrap
#' @export
resid_bootstrap.lme <- function(model, .f, B){
  
  .f <- match.fun(.f)
  
  setup <- .setup.lme(model, type = "residual")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.resids.lme(
        b = setup$b, 
        e = setup$e, 
        Xbeta = setup$Xbeta,
        Zlist = setup$Zlist
      )
    )
  )
  
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(updated.model(model = model, new.y = x))
  })
  
  
  lmeresampler:::.bootstrap.completion(model, tstar, B, .f, type = "residual")
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
  
  y.star <- as.data.frame(y.star)
  
  # Extract bootstrap statistics
  if(reb_type == 2) .f <- extract_parameters.lme
  
  tstar <- purrr::map_dfc(y.star, function(x) {
    .f(updated.model(model = model, new.y = x))
  })
  
  # Extract original statistics
  t0 <- .f(model)
  
  # Postprocessing for REB/2
  if(reb_type == 2) 
    tstar <- lmeresampler:::.postprocess.reb2(t0, tstar, nbeta = length(fixef(model)), B = B)
  
  # Format for return
  lmeresampler:::.bootstrap.completion(model, tstar, B, .f, type = paste("reb", reb_type, sep = ""))
}





#' @rdname cgr_bootstrap
#' @inheritParams bootstrap
#' @export
cgr_bootstrap.lme <- function(model, .f, B){
  
  .f <- match.fun(.f)
  
  setup <- .setup.lme(model, type = "cgr")
  
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
  
  tstar <- purrr::map_dfc(ystar, function(x) {
    .f(updated.model(model = model, new.y = x))
  })
  
  
  lmeresampler:::.bootstrap.completion(model, tstar, B, .f, type = "cgr")
}



