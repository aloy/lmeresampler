#' @title Print a summary of an \code{lmeresamp} object
#'
#' @description
#' Print summary statistics and confidence intervals, if desired, for an \code{lmeresamp} object.
#'
#' @details
#' If the bootstrap statistics are stored in a vector (as opposed to a data frame or tibble), 
#' then summary statistics will be calculated and printed. The printed data frame will include
#' the name of the term (if applicable), the observed value (\code{observed}), the mean of the bootstrap replicated 
#' (\code{rep.mean}), the standard error (\code{se}), and the bootstrap bias estimate (\code{bias}).
#' In addition, the number of resamples will be printed. If any messages, warnings, or errors were
#' generated during the bootstrap procedure, they will be summarized below, and you should check the 
#' \code{message}, \code{warning}, and \code{error} elements of the \code{lmeresamp} object to
#' investigate further.
#'
#' @param x The lmeresamp object to print.
#' @param ci A logical value specifying whether confidence intervals should be printed.
#' @param ... not used
#'
#' @rdname print
#' @export 
#' @method print lmeresamp
print.lmeresamp <- function(x, ci = FALSE, ...){
  
  summary.lmeresamp(x)
  
  if(ci == TRUE){
    cat(paste("\n"))
    cat(paste("\n"))
    confint(x)
  }
}



#' @title Print coefficients from a \code{coef_tbl} object
#' 
#' @description
#' Print table of coefficients produced by \code{bootstrap_pvals}.
#' 
#' @param x the coef_tbl object to print
#' @method print coef_tbl
#' @param ... not used
#' @export 
print.coef_tbl <- function(x, ...) {
  cat(paste("Bootstrap type:", x$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", x$B, "\n"))
  cat(paste("\n"))
  print(x$coefficients)
  
}