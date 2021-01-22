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