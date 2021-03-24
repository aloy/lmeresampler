#' @title summary
#'
#' @description
#' Summarizes the output of the bootstrap call.
#'
#' @details
#' This function is given \code{object, ci} and uses them to print the bootstrap result
#' and confidence intervals if the user wishes to have them.
#'
#' @param object ci A logical value specifying whether confidence intervals should be printed.
#' @param ... not used
#'
#' @rdname summary
#' @export 
#' @method summary lmeresamp
summary.lmeresamp <- function(object, ci = FALSE, ...){
  
  cat(paste("Bootstrap type:", object$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", object$R, "\n"))
  cat(paste("\n"))
  print(object$stats)
  
  if(ci == TRUE){
    confint(object)
  }
}
