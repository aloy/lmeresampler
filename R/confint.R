#' @title confint
#'
#' @description
#' Performs confidence intervals on an lmeresamp object.
#'
#' @details
#' This function is given \code{object, method, level} and uses them to calculate the
#' confidence interval(s).
#'
#' @param object The lmeresamp object for which confidence intervals should be computed.
#' @param method The type of confidence intervals that should be executed.
#' @param level The level at which the confidence interval should be calculated.
#'
#' @rdname confint
#' @export 
# bootstrap CI method for object of class lmeresamp
confint.lmeresamp <- function(object, method, level) {
  library(dplyr)
  
  if(missing(level)){
    level <- 0.95
  } else if(!level > 0 && !level < 1){
    stop("please specify a confidence level between 0 and 1")
  }
  
  if(missing(method)){
    method <- "all"
  } else if(!method %in% c("norm", "boot-t", "perc", "all")) {
    stop("'method' must be either 'norm', 'boot-t', 'perc', or 'all'")
  }
  
  if(class(object$model) == "lmerMod"){ ## normal t
    
    if(method == "norm") {
      
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      colnames(con)[1] <- "norm.t.lower"
      norm.t.lower <- con[1]
      colnames(con)[2] <- "norm.t.upper"
      norm.t.upper <- con[2]
      
      .norm.t.completion(norm.t.lower, norm.t.upper, level)
      
    } else if(method == "boot-t"){ ## boot t
      ## note, this is for fixef only, regardless of what .f is
      
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- lme4::getME(object$model, "beta")
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$coefficients[, 2] # fixef
      
      ### table of estimates and sds for boot_t calculation
      t.stats <- cbind(model.fits, model.sds)
      
      .boot.t.completion(object, level, t.stats)
      
    } else if(method == "perc"){ ## percentile t
      
      .perc.t.completion(object, level)
      
    } else if(method == "all"){
      
      ## normal-t
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      colnames(con)[1] <- "norm.t.lower"
      norm.t.lower <- con[1]
      colnames(con)[2] <- "norm.t.upper"
      norm.t.upper <- con[2]
      
      cat(paste("\n"))
      .norm.t.completion(norm.t.lower, norm.t.upper, level)
      
      ## boot-t
      ### this gets all of the means (estimates) for boot_t calculation
      con <- data.frame(lme4::confint.merMod(object$model, level = level))
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- lme4::getME(object$model, "beta")
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$coefficients[, 2] # fixef
      
      ### table of estimates and sds for boot_t calculation
      t.stats <- cbind(model.fits, model.sds)
      
      .boot.t.completion(object, level, t.stats)
      
      ## percentile-t
      .perc.t.completion(object, level)
      
    } 
  } else if(class(object$model) == "lme"){
    
    if(method == "norm"){ ## normal t
      
      ### intervals returns the same estimates as doing summary() so let's use it!
      con <- nlme::intervals(object$model, level = level, which = "all")
      norm.t.lower <- con$fixed[, 1]  # lower interval for all fixedef
      norm.t.lower <- append(norm.t.lower, con$reStruct[[1]][1]) # lower interval for ranef
      norm.t.lower <- append(norm.t.lower, con$sigma[1]) # lower interval for sigma
      norm.t.lower <- unlist(norm.t.lower)
      names(norm.t.lower)[(length(con$fixed[, 1]) + 1) : (length(con$fixed[, 1]) + length(con$reStruct[[1]][1]))] <- row.names(con$reStruct[[1]][1])
      names(norm.t.lower)[(length(con$fixed[, 1]) + length(con$reStruct[[1]][1]) + 1)] <- "sigma"
      
      norm.t.upper <- con$fixed[, 3]  # upper interval for all fixedef
      norm.t.upper <- append(norm.t.upper, con$reStruct[[1]][3]) # upper interval for ranef
      norm.t.upper <- append(norm.t.upper, con$sigma[3]) # upper interval for sigma
      norm.t.upper <- unlist(norm.t.upper)
      names(norm.t.upper)[(length(con$fixed[, 3]) + 1) : (length(con$fixed[, 3]) + length(con$reStruct[[1]][3]))] <- row.names(con$reStruct[[1]][3])
      names(norm.t.upper)[(length(con$fixed[, 3]) + length(con$reStruct[[1]][3]) + 1)] <- "sigma"
      
      .norm.t.completion(norm.t.lower, norm.t.upper, level)
      
    } else if(method == "boot-t"){  ## boot t
      con <- nlme::intervals(object$model, level = level, which = "all")
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- con$fixed[, 2]  # estimates for all fixef
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$tTable[, 2] # fixef 

      t.stats <- cbind(model.fits, model.sds)
      
      .boot.t.completion(object, level, t.stats)
      
    } else if(method == "perc"){ ## percentile t
      
      .perc.t.completion(object, level)
      
    } else if(method == "all"){
      
      ## normal-t
      ### intervals returns the same estimates as doing summary() so let's use it!
      con <- nlme::intervals(object$model, level = level, which = "all")
      norm.t.lower <- con$fixed[, 1]  # lower interval for all fixedef
      norm.t.lower <- append(norm.t.lower, con$reStruct[[1]][1]) # lower interval for ranef
      norm.t.lower <- append(norm.t.lower, con$sigma[1]) # lower interval for sigma
      norm.t.lower <- unlist(norm.t.lower)
      names(norm.t.lower)[(length(con$fixed[, 1]) + 1) : (length(con$fixed[, 1]) + length(con$reStruct[[1]][1]))] <- row.names(con$reStruct[[1]][1])
      names(norm.t.lower)[(length(con$fixed[, 1]) + length(con$reStruct[[1]][1]) + 1)] <- "sigma"
      
      norm.t.upper <- con$fixed[, 3]  # upper interval for all fixedef
      norm.t.upper <- append(norm.t.upper, con$reStruct[[1]][3]) # upper interval for ranef
      norm.t.upper <- append(norm.t.upper, con$sigma[3]) # upper interval for sigma
      norm.t.upper <- unlist(norm.t.upper)
      names(norm.t.upper)[(length(con$fixed[, 3]) + 1) : (length(con$fixed[, 3]) + length(con$reStruct[[1]][3]))] <- row.names(con$reStruct[[1]][3])
      names(norm.t.upper)[(length(con$fixed[, 3]) + length(con$reStruct[[1]][3]) + 1)] <- "sigma"
      
      cat(paste("\n"))
      .norm.t.completion(norm.t.lower, norm.t.upper, level)
      
      ## boot-t
      con <- nlme::intervals(object$model, level = level, which = "all")
      
      ### this gets all of the means (estimates) for boot_t calculation
      model.fits <- con$fixed[, 2]  # estimates for all fixef
      
      ### sd of of estimates for boot_t calculation
      out <- summary(object$model)
      model.sds <- out$tTable[, 2] # fixef 
      
      t.stats <- cbind(model.fits, model.sds)
      
      .boot.t.completion(object, level, t.stats)
      
      ## percentile t
      .perc.t.completion(object, level)
    }
  }
}


#' @title Percentile-t interval completion
#'
#' @description
#' Execute the percentile-t interval process
#'
#' @details
#' This function uses \code{object} and \code{level} to calculate a percentile-t
#' interval for the fixed and random components
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.perc.t.completion <- function(object, level){
  
  perc.t.lower <- apply(object$replicates, 2, function(x) {
    round(quantile(x, (1 - level)/2), 8)
  })
  
  perc.t.upper <- apply(object$replicates, 2, function(x) {
    round(quantile(x, level + (1 - level)/2), 8)
  })
  
  perc.t <- data.frame(cbind(perc.t.lower, perc.t.upper))
  rownames(perc.t) <- rownames(object$stats)
  
  conf.lev <- level*100
  cat(paste(conf.lev, "percent percentile-t interval: \n"))
  print(perc.t)
  cat(paste("\n"))
}

#' @title Bootstrap-t interval completion
#'
#' @description
#' Finish the bootstrap-t interval process
#'
#' @details
#' This function uses \code{object} and \code{level} to calculate a bootstrap-t
#' interval for the fixed and random components
#'
#' @param object An lmeresamp object
#' @param level A confidence level
#'
#' @keywords internal
#' @noRd
.boot.t.completion <- function(object, level, t.stats){
  
  numFixef <- length(rownames(t.stats))
  
  t.stats <- as.data.frame(t.stats)
  t.vals <- t.stats %>% 
    mutate(boot_t  = (object$stats$rep.mean[1:numFixef] - t.stats[, 1])/(object$stats$se[1:numFixef]/sqrt(object$R))) %>%
    mutate(boot.t.lower = (t.stats[, 1] - quantile(boot_t, level + (1 - level)/2) * t.stats[, 2]/sqrt(object$R))) %>%
    mutate(boot.t.upper = (t.stats[, 1] - quantile(boot_t, (1 - level)/2) * t.stats[, 2]/sqrt(object$R)))
  
  boot.t <- t.vals %>%
    select(boot.t.lower, boot.t.upper)
  rownames(boot.t) <- rownames(t.stats)
  
  conf.lev <- level*100
  cat(paste(conf.lev, "percent bootstrap-t interval: \n"))
  print(boot.t)
  cat(paste("\n"))
}

#' @title Normal-t interval completion
#'
#' @description
#' Complete the normal-t interval process
#'
#' @details
#' This function uses \code{norm.t.lower} and \code{norm.t.upper} to make a normal-t
#' interval for the fixed and random components
#'
#' @param norm.t.lower The lower bound of the interval
#' @param norm.t.upper The upper bound of the interval
#'
#' @keywords internal
#' @noRd
.norm.t.completion <- function(norm.t.lower, norm.t.upper, level){
  
  norm.t <-  data.frame(cbind(norm.t.lower, norm.t.upper))
  
  conf.lev <- level*100
  cat(paste(conf.lev, "percent normal-t interval: \n"))
  print(norm.t)
  cat(paste("\n"))
}

