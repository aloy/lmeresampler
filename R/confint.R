#' @title Calculate confidence intervals for a lmeresamp object
#'
#' @description
#' Calculates normal, basic, and percentile bootstrap confidence intervals 
#' from a \code{lmeresamp} object.
#'
#'
#' @param object The lmeresamp object for which confidence intervals should be computed.
#' @param type A character string giving the type of confidence intervals that should be calculated. 
#' This should be a subset of \code{c("norm", "basic", "perc")} (for normal, basic, and percentile
#' bootstrap confidence intervals, respectively), or \code{"all"}.
#' @param level The level at which the confidence interval should be calculated. 
#' @param parm not used
#' @param ... not used
#' 
#' @return 
#' A tibble with columns term, estimate, lower, upper, type, and level.
#' 
#'
#' @rdname confint
#' @export 
#' @importFrom stats qnorm
confint.lmeresamp <- function(object, parm, level = 0.95, 
                              type = c("all", "norm", "basic", "perc"), 
                              ...) {
  term <- estimate <- lower <- upper <- NULL
  
  if(!level > 0 && !level < 1){
    stop("please specify a confidence level between 0 and 1")
  }
  
  type <- match.arg(type)
  
  terms <- names(object$observed)
  if(is.null(terms)) terms <- ""
  orig <- dplyr::tibble(term = terms, estimate = object$observed)

  ci.out <- NULL
  if(any(type == "all" | type == "norm")) {
    ci.out <- c(ci.out, list(.norm_ci(object, level)))
  }
  
  if(any(type == "all" | type == "basic")) {
    ci.out <- c(ci.out, list(.basic_ci(object, level)))
  }
  
  if(any(type == "all" | type == "perc")) {
    ci.out <- c(ci.out, list(.perc_ci(object, level)))
  }
  
  ci.out <- lapply(
    ci.out, 
    function(x) dplyr::bind_cols(orig, lower = x[,1], upper = x[,2])
  )
  
  if(type == "all") type <- c("norm", "basic", "perc")
  names(ci.out) <- type
  
  dplyr::bind_rows(ci.out, .id = "type") %>% 
    dplyr::mutate(level = level) %>% 
    dplyr::select(term, estimate, lower, upper, dplyr::everything())
}



#' Calculate Percentile Bootstrap CI
#' 
#' @description
#' Calculate a percentile bootstrap interval
#' 
#' @details
#' This function uses \code{object} and \code{level} to calculate a percentile
#' interval for the components specified by .f
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.perc_ci <- function(object, level){
  if(typeof(object$replicates) == "list") {
    t(
      apply(object$replicates, 2, quantile, probs = (1 + c(-level, level)) / 2, na.rm = TRUE)
    )
  } else {
    t(quantile(object$replicates, probs = (1 + c(-level, level)) / 2, na.rm = TRUE))
  }
}


#' Calculate Basic Bootstrap CI
#' 
#' @description
#' Calculate a basic bootstrap interval
#' 
#' @details
#' This function uses \code{object} and \code{level} to calculate a basic
#' interval for the components specified by .f
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.basic_ci <- function(object, level){
  if(typeof(object$replicates) == "list") {
    quants <- apply(object$replicates, 2, quantile, probs = (1 + c(level, -level))/2)
    ci <- 2 * object$observed - t(quants)
    colnames(ci) <- rev(colnames(ci))
  } else{
    quants <- quantile(object$replicates, probs = (1 + c(level, -level))/2)
    ci <- 2 * object$observed - t(quants)
  }
  ci
}


#' Calculate Normal Bootstrap CI
#' 
#' @description
#' Calculate a normal bootstrap interval
#' 
#' @details
#' This function uses \code{object} and \code{level} to calculate a normal
#' interval for the components specified by .f
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.norm_ci <- function(object, level){
  if(typeof(object$replicates) == "list") {
    se   <- apply(object$replicates, 2, sd, na.rm = TRUE)
  } else {
    se <- sd(object$replicates, na.rm = TRUE)
  }
  merr <- se * qnorm((1 + level) / 2)
  bias <- object$stats$bias
  orig <- object$observed
  ci <- cbind(orig - bias - merr, orig - bias + merr)
  colnames(ci) <- paste0(100 * (1 + c(-level, level))/2, "%")
  ci
}
