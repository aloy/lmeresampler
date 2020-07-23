#' @rdname print
#' @export
print.lmeresamp <- function(x, ci = FALSE){
  
  if(x$type = "reb"){
    cat(paste("Bootstrap type: REB", x$reb_type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", x$B, "\n"))
    cat(paste("\n"))
    print(x$stats)
    
    if(ci == TRUE){
      confint.lmeresamp(x)
    }
  }
  else{
    cat(paste("Bootstrap type:", x$type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", x$B, "\n"))
    cat(paste("\n"))
    print(x$stats)
    
    if(ci == TRUE){
      confint.lmeresamp(x)
    }
  }
}
