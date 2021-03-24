#' @title summary
#'
#' @description
#' Summarizes the output of the bootstrap call.
#'
#' @details
#' This function is given \code{object, ci} and uses them to print the bootstrap result
#' and confidence intervals if the user wishes to have them.
#'
#' @param object The lmeresamp object to be summarized.
#' @param ... not used
#'
#' @rdname summary
#' @export 
#' @method summary lmeresamp
summary.lmeresamp <- function(object, ...){
  
  # use this to test warnings
  # lme_res_boot <- bootstrap(vcmodB, .f = fixef, type = "residual", B = 100)
  
  cat(paste("Bootstrap type:", object$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", object$R, "\n"))
  cat(paste("\n"))
  print(object$stats)
  cat(paste("\n"))
  
  message_count <- 0
  for(i in length(object$message)){
    if(!is.null(object$message[i])){
      message_count <- message_count + 1
    }
  }
  
  warning_count <- 0
  for(j in length(object$warning)){
    if(!is.null(object$warning[j])){
      warning_count <- warning_count + 1
    }
  }
  
  error_count <- 0
  for(k in length(object$error)){
    if(!is.null(object$error[k])){
      error_count <- error_count + 1
    }
  }
  
  cat(paste("There were", message_count, "messages,", warning_count, "warnings, and", error_count, "errors."))
  
}
