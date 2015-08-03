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
#' @param replace
#' @param reb_type Chooses the type of REB bootstrap
#' 
#' @references
#'    Carpenter:2003uy
#'    Chambers:2013ba
#'    Morris:2002tj
#'    vanderLeeden:208kv
bootstrap <- function(model, fn, type, B, replace, reb_type) {
  UseMethod("bootstrap", model)
}


#' @title Parametric Bootstrap
#'
#' @description
#' The Parametric Bootstrap is uses the parametrically estimated
#' distribution function of the data to generate bootstrap samples.
#'
#' @details
#' This function extracts the fixed effects, simulates from the model, refits the model
#' and then returns the results in a list.
#'
#' @export
#' @inheritParams bootstrap
#'
#' @return list
#'
#' @references
#' Chambers:2013ba
#' vanderLeeden:208kv
parametric_bootstrap <- function(model, fn, B) {
  UseMethod("parametric_bootstrap", model)
}

#' @title Residual Bootstrap
#'
#' @description
#' The Residual Bootstrap uses residuals to generate bootstrap samples.
#'
#' @details
#' Resamples the residuals and complete the bootstrap process.
#'
#' @export
#' @inheritParams bootstrap
#'
#' @return list
#'
#' @references
#' vanderLeeden:208kv
resid_bootstrap <- function(model, fn, B) {
  UseMethod("resid_bootstrap", model)
}

#' @title Cases Bootstrap
#'
#' @description
#' The Cases Bootstrap samples entire cases to generate the bootstrap.
#'
#' @details
#' add details
#' 
#' @export
#' @inheritParams bootstrap
#'
#' @return list
#'
#' @references
#' vanderLeeden:208kv
case_bootstrap <- function(model, fn, B, replace) {
  UseMethod("case_bootstrap", model)
}

#' CGR Bootstrap
#'
#' @description
#' CGR Bootstrap
#'
#' @details
#' add details later
#'
#' @inheritParams bootstrap
#'
#' @return list
#'
#' @references Carpenter, J. R., Goldstein, H., and Rasbash, J. (2003)
#'  A novel bootstrap procedure for assessing the relationship 
#'  between class size and achievement. \emph{Journal of the Royal 
#'  Statistical Society. Series C. Applied Statistics}, 52(4), 431--443. 
#'  doi:10.1111/1467-9876.00415
cgr_bootstrap <- function(model, fn, B) {
  UseMethod("cgr_bootstrap", model)
}

#' @title REB Bootstrap
#'
#' @description
#' REB Bootstrap
#'
#' @details
#' add details
#'
#' @export
#' @inheritParams bootstrap
#'
#' @return list
#'
#' @references
#' Chambers, R. and Chandra, H. (2013) 
#' A Random Effect Block Bootstrap for Clustered Data. 
#' \emph{Journal of Computational and Graphical Statistics}, 
#' 22(2), 452–470. doi:10.1080/10618600.2012.681216
reb_bootstrap <- function(model, fn, B, reb_type) {
  UseMethod("reb_bootstrap", model)
}