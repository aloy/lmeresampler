#' @title Calculate summary statistics for \code{lmeresamp} objects
#'
#' @description
#' Calculate basic summary statistics such as the mean, standard error, and bias of the
#' bootstrap replicates.
#'
#' @details
#' If the bootstrap statistics are stored in a vector (as opposed to a data frame or tibble), 
#' then summary statistics will be calculated. The printed data frame will include
#' the name of the term (if applicable), the observed value (\code{observed}), the mean of the bootstrap replicated 
#' (\code{rep.mean}), the standard error (\code{se}), and the bootstrap bias estimate (\code{bias}).
#' In addition, the number of resamples will be printed. If any messages, warnings, or errors were
#' generated during the bootstrap procedure, they will be summarized below, and you should check the 
#' \code{message}, \code{warning}, and \code{error} elements of the \code{lmeresamp} object to
#' investigate further.
#'
#' @param object The lmeresamp object to be summarized.
#' @param ... not used
#'
#' @rdname summary
#' @export 
#' @method summary lmeresamp
summary.lmeresamp <- function(object, ...){
  
  cat(paste("Bootstrap type:", object$type, "\n"))
  cat(paste("\n"))
  cat(paste("Number of resamples:", object$B, "\n"))
  cat(paste("\n"))
  print(as.data.frame(object$stats))
  
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
  
  top_message <- names(sort(summary(object$message), decreasing=T)[1])
  if(!is.null(top_message)){
    cat(paste("\n"))
    cat(paste("The most commonly occurring message was:", top_message))
    cat(paste("\n"))
  }
  
  top_warning <- names(sort(summary(object$warning), decreasing=T)[1])
  if(!is.null(top_warning)){
    cat(paste("\n"))
    cat(paste("The most commonly occurring warning was:", top_warning))
    cat(paste("\n"))
  }
  
  top_error <- names(sort(summary(object$error), decreasing=T)[1])
  if(!is.null(top_error)){
    cat(paste("\n"))
    cat(paste("The most commonly occurring error was:", top_error))
  }
  
}
