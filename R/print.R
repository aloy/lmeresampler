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
#'
#' @rdname print
#' @export 
#' @method print lmeresamp
print.lmeresamp <- function(x, ci = FALSE){
  
  if(is.null(x$type)){
    cat(paste("Bootstrap type: REB", x$reb_type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", x$R, "\n"))
    cat(paste("\n"))
    print(x$stats)
    
    if(ci == TRUE){
      confint(x)
    }
  } else{
    cat(paste("Bootstrap type:", x$type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", x$R, "\n"))
    cat(paste("\n"))
    print(x$stats)
    
    if(ci == TRUE){
      confint(x)
    }
  }
}
