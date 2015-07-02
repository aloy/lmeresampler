#' @rdname bootstrap
bootstrap.lme <- function (model, fn, type, B, extra_step, reb_type){
  switch(type,
         par = parametric_bootstrap.lme(model, fn, B),
         resid = resid_bootstrap.lme(model, fn, B),
         case = case_bootstrap.lme(model, fn, B, extra_step = FALSE),
         cgr = cgr_bootstrap.lme(model, fn, B),
         reb = reb_bootstrap.lme(model, fn, B, reb_type = 0))
}


#' @rdname parametric_bootstrap
#' @export
parametric_bootstrap.lme <- function(model, fn, B){
  # getVarCov.lme is the limiting factor...
  if (length(model$group) > 1) 
    stop("not implemented for multiple levels of nesting")
  
  # Match function
  fn <- match.fun(fn)
  
  # Extract fixed effects
  model.fixef <- nlme::fixef(model)
  
  ystar <- nlmeU::simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  ystar <- data.frame(ystar)
  
  t0 <- fn(model)
  
  # t.res <- list()
  # length(t.res) <- B
  #    t.res <- matrix(0, ncol = 2, nrow = B)
  #   for(i in 1:B){
  #     myin <- ystar[,i]
  #     model.update <- update(object = model, fixed = myin ~ .)
  #     t.res[i,] <- fn(model.update)
  #   }
  
  
  
  res <- lapply(ystar, function(y) {
    fit <- tryCatch(fn(updated.model(model = model, new.y = y)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
#   for(i in 1:B){
#     #     myin <- ystar[,i]
#     #     model.update <- nlme:::update.lme(object = model, fixed = myin ~ .)
#     #     t.res[i,] <- fn(model.update)
#     fit <- tryCatch(fn(updated.model(model = model, new.y = ystar[,i])),  error = function(e) e)
#     if (inherits(fit, "error")) {
#       structure(rep(NA, length(t0)), fail.msgs = fit$message)
#     } else{
#       t.res[[i]] <- fit
#     }
    ### NOTE: Need to check if there was an error in try, and then use fn()
    # tmp.mod <- updated.model(model = model, new.y = ystar[,i])
    # t.res[[i]] <- fn(try.fit)
  # }
  tstar <- do.call('cbind', res)
  # tstar <- data.frame(tstar)
  
  #   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
  #   
  #   tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}



#' @rdname case_bootstrap
#' @export
case_bootstrap.lme <- function (model, fn, B, replace){
  
  data <- model$data
  data$.id <- seq_len(nrow(data))
  clusters <- c(names(model$groups), ".id")
  
  ## ADD ERROR CHECKS!!
  
  t0 <- fn(model)
  
  rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(dat = data, cluster = clusters, replace = replace))))
  
  res <- lapply(rep.data, function(df) {
    fit <- tryCatch(fn(updated.model(model = model, new.data = df)),  
                    error = function(e) e)
    if (inherits(fit, "error")) {
      structure(rep(NA, length(t0)), fail.msgs = fit$message)
    } else{
      fit
    }
  }
  )
  
  tstar <- do.call('cbind', res)
  
  rownames(tstar) <- names(t0)
  colnames(tstar) <- names(res) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  
  if ((numFail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", B, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  } else fail.msgs <- NULL
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  attr(RES, "bootFail") <- numFail
  attr(RES, "boot.fail.msgs") <- fail.msgs
  return(RES)
}