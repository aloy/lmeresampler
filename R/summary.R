#' @title summary
#'
#' @description
#' Summarizes the output of the bootstrap call.
#'
#' @details
#' This function is given \code{x, ci} and uses them to print the bootstrap result
#' and confidence intervals if the user wishes to have them.
#'
#' @param x The lmeresamp object to print.
#' @param ci A logical value specifying whether confidence intervals should be printed.
#' @param ... not used
#'
#' @rdname summary
#' @export 
#' @method summary lmeresamp
summary.lmeresamp <- function(x, ci = FALSE, ...){
  
  cat(paste("Bootstrap type:", x$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", x$R, "\n"))
  cat(paste("\n"))
  print(x$stats)
  
  if(ci == TRUE){
    confint(x)
  }
}
