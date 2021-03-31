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
  
  cat(paste("Bootstrap type:", x$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", x$B, "\n"))
  cat(paste("\n"))
  print(as.data.frame(x$stats))
  
  message_count <- 0
  for(i in length(x$message)){
    if(!is.null(x$message[[i]])){
      message_count <- message_count + 1
    }
  }
  
  warning_count <- 0
  for(j in length(x$warning)){
    if(!is.null(x$warning[[j]])){
      warning_count <- warning_count + 1
    }
  }
  
  error_count <- 0
  for(k in length(x$error)){
    if(!is.null(x$error[[k]])){
      error_count <- error_count + 1
    }
  }
  
  cat(paste("\n"))
  cat(paste("There were", message_count, "messages,", warning_count, "warnings, and", error_count, "errors."))
  cat(paste("\n"))
  
  if(ci == TRUE){
    confint(x)
  }
}
