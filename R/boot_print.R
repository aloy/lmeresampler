# print function that will be called in bootstrap_completion()
boot_print(type = type, B = B, ci = FALSE) <- function(x){
  
  if(type = "reb"){
    cat(paste("Bootstrap type: REB", reb_type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", B, "\n"))
    cat(paste("\n"))
    print(stats)
    
    # Confidence intervals from ci_prep
    # put them all together 
    norm.t.cis <-  data.frame(cbind(norm.t.lower, norm.t.upper))
    other.cis <- data.frame(cbind(boot.t, perc.t.lower, perc.t.upper))
    cat(paste("\n"))
    cat(paste("95% normal t-interval: \n"))
    print(norm.t.cis)
    cat(paste("\n"))
    cat(paste("95% bootstrap-t and percentile confidence intervals: \n"))
    print(other.cis)
    
  }
  else{
    cat(paste("Bootstrap type:", type, "\n"))
    cat(paste("\n"))
    cat(paste("Number of resamples:", B, "\n"))
    cat(paste("\n"))
    print(stats)
    
    # Confidence intervals from ci_prep
    # put them all together 
    norm.t.cis <-  data.frame(cbind(norm.t.lower, norm.t.upper))
    other.cis <- data.frame(cbind(boot.t, perc.t.lower, perc.t.upper))
    cat(paste("\n"))
    cat(paste("95% normal t-interval: \n"))
    print(norm.t.cis)
    cat(paste("\n"))
    cat(paste("95% bootstrap-t and percentile confidence intervals: \n"))
    print(other.cis)
  }
}