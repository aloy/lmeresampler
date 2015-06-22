#' @title Bootstrap Nested LMEs
#'
#' @description
#' \code{bootstrap} helps streamline the bootstrap process for the parametric,
#' residual, cases, CGR, and REB bootstraps.
#'
#'
#' @export
#' @param model The model object you wish to bootstrap.
#' @param fn A function returning the statistic(s) of interest.
#' @param type A character string indicating the type of bootstrap that is being
#'    requested. Possible values are \code{"par"} (parametric), \code{"resid"} 
#'    (residual), \code{"case"}, \code{"cgr"} (), or one of the three versions 
#'    of the random effect block bootstrap ("reb", "reb1", "reb2").
#' @param B The number of bootstrap resamples.
#' 
#' @references
#'    Carpenter:2003uy
#'    Chambers:2013ba
#'    Morris:2002tj
#'    vanderLeeden:208kv
bootstrap <- function(model, fn, type, B) {
  UseMethod("boostrap", model)
}