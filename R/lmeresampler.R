#' lmeresampler: A package for bootstrapping nested linear mixed-effects models
#'
#' The \pkg{lme4} and \pkg{nlme} packages have made fitting nested 
#' linear mixed-effects (LME) models quite easy. Using the the 
#' functionality of these packages we can easily use maximum 
#' likelihood or restricted maximum likelihood to fit a 
#' model and conduct inference using our parametric toolkit. 
#' In practice, the assumptions of our model are often violated 
#' to such a degree that leads to biased estimators and 
#' incorrect standard errors. In these situations, resampling 
#' methods such as the bootstrap can be used to obtain consistent 
#' estimators and standard errors for inference. 
#' \code{lmeresampler} provides an easy way to bootstrap nested 
#' linear-mixed effects models using either fit using either \pkg{lme4} or 
#' \pkg{nlme}.
#' 
#' 
#' A variety of bootstrap procedures are available: 
#' \itemize{
#'    \item the parametric bootstrap: \code{\link{parametric_bootstrap}}
#'    \item the residual bootstrap: \code{\link{resid_bootstrap}}
#'    \item the cases (i.e. non-parametric) bootstrap: \code{\link{case_bootstrap}}
#'    \item the random effects block (REB) bootstrap: \code{\link{reb_bootstrap}}
#'    \item the Wild bootstrap: \code{\link{wild_bootstrap}}
#' }
#' 
#' In addition to the individual bootstrap functions, \code{lmeresampler} provides
#' a unified interface to bootstrapping LME models in its  \code{bootstrap} function.
#' @docType package
#' @name lmeresampler
#' @aliases lmeresampler package-lmeresampler
#' @keywords package
NULL