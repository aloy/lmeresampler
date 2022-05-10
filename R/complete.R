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
.bootstrap.completion <- function(model, tstar, B, .f, type = type, warnings){
  
  # tstar <- do.call("cbind", tstar) # Can these be nested?
  # row.names(tstar) <- names(t0)
  
  # if((nfail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
  #   warning("some bootstrap runs failed (", nfail, "/", nsim, ")")
  #   fail.msgs <- purrr::map_chr(tstar[bad.runs], .f = attr, FUN.VALUE = character(1),
  #                               "fail.msgs")
  # } else fail.msgs <- NULL
  
  # prep for stats df

  t0 <- .f(model)
  t0len <- length(t0)
  nsim <- length(tstar)
  observed <- t0
  
  if(is.numeric(t0)) {
    
    replicates <- dplyr::bind_rows(tstar)
    rep.mean <- colMeans(replicates[, 1:t0len])
    se <- unlist(purrr::map(replicates[, 1:t0len], sd))
    bias <- rep.mean - observed
    
    if(t0len == 1) {
      stats <- dplyr::tibble(observed, rep.mean, se, bias)
    } else{
      # Check for names
      nms <- unlist(lapply(tstar, names))
      if(is.null(nms)) 
        warning("Lists of unnamed vectors are converted to data frames.\nPlease create named vectors in .f() if this is not the desired behavior.",
                call. = FALSE)
      stats <- dplyr::tibble(term = names(t0), observed, rep.mean, se, bias)
    }
    
  } else{
    if(is.data.frame(t0)) {
      .ids <- rep(seq_along(tstar), times = vapply(tstar, nrow, FUN.VALUE = 0L))
      replicates <- dplyr::bind_rows(tstar) %>% dplyr::mutate(.n = .ids)
    }
    stats <- NULL
  }

  
  if (inherits(model, "lme")) data <- model$data
  else data <- model@frame
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, B = B, data = data,
                        seed = .Random.seed, type = type, call = match.call(),
                        message = warnings$message, warning = warnings$warning, error = warnings$error), 
                   class = "lmeresamp")
  
  # attr(RES,"bootFail") <- nfail
  # attr(RES,"boot.fail.msgs") <- fail.msgs
  return(RES)
}
