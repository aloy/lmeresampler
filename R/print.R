#' @title print
#'
#' @description
#' Prints the output of the bootstrap call.
#'
#' @details
#' This function is given \code{x, ci} and uses them to print the bootstrap result
#' and confidence intervals if the user wishes to have them.
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
    confint(x)
  }
}
