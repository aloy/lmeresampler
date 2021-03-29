#' @title summary
#'
#' @description
#' Summarizes the output of the bootstrap call.
#'
#' @details
#' This function is given \code{object} and uses them to summarize the bootstrap result.
#'
#' @param object The lmeresamp object to be summarized.
#' @param ... not used
#'
#' @rdname summary
#' @export 
#' @method summary lmeresamp
summary.lmeresamp <- function(object, ...){
  
  print.lmeresamp(object)
  
  message_count <- 0
  for(i in length(object$message)){
    if(!is.null(object$message[[i]])){
      message_count <- message_count + 1
    }
  }
  
  warning_count <- 0
  for(j in length(object$warning)){
    if(!is.null(object$warning[[j]])){
      warning_count <- warning_count + 1
    }
  }
  
  error_count <- 0
  for(k in length(object$error)){
    if(!is.null(object$error[[k]])){
      error_count <- error_count + 1
    }
  }
  
  cat(paste("\n"))
  cat(paste("There were", message_count, "messages,", warning_count, "warnings, and", error_count, "errors."))
  cat(paste("\n"))
  
  # finding most commonly occuring message/warning/error
  object$message <- as.factor(object$message)
  object$warning <- as.factor(object$warning)
  object$error <- as.factor(object$error)
  
  top_message <- names(sort(summary(object$message), decreasing=T)[[1]])
  if(!is.null(top_message)){
    cat(paste("The most commonly occuring message was:"))
    cat(paste("\n"))
  }
  
  top_warning <- names(sort(summary(object$warning), decreasing=T)[[1]])
  if(!is.null(top_warning)){
    cat(paste("The most commonly occuring warning was:"))
    cat(paste("\n"))
  }
  
  top_error <- names(sort(summary(object$error), decreasing=T)[[1]])
  if(!is.null(top_error)){
    cat(paste("The most commonly occuring error was:"))
  }
  
}
