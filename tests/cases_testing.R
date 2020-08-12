library(tidyverse)
library(lme4)
# load in data
deer <- read.csv("http://math.carleton.edu/Chihara/Stat330/Deer.csv")

# clean up data a bit
library(dplyr)
deer <- deer %>% mutate(prop = PositiveCervi/Sampled)
deer <- deer %>% mutate(Id = row_number())

deer.glmer8 <- glmer(prop ~ Open + Scrub + I(Scrub^2) + (1|Id) , data = deer,
                     family = binomial, weights = Sampled)

model <- deer.glmer8

#' @rdname case_bootstrap
#' @export
case_bootstrap.lmerMod <- function(model, .f, B, resample, type){
  
  data <- model@frame
  # data$.id <- seq_len(nrow(data))
  clusters <- c(rev(names(lme4::getME(model, "flist"))), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping.")
  
  # rep.data <- purrr::map(integer(B), function(x) .cases.resamp(model = model, dat = data, cluster = clusters, resample = resample))
  tstar <- purrr::map(integer(B), function(x) .cases.resamp(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  RES <- .bootstrap.completion(model, tstar, B, .f, type)
  return(RES)
}

# # Using recursion allows for a very general function...
# # How can I speed this up?
# .cases.resamp <- function(dat, cluster, resample) {
#   # exit early for trivial data
#   if(nrow(dat) == 1 || all(resample==FALSE))
#     return(dat)
#   
#   # sample the clustering factor
#   cls <- sample(unique(dat[[cluster[1]]]), resample=resample[1])
#   
#   # subset on the sampled clustering factors
#   sub <- lapply(cls, function(b) dat[dat[[cluster[1]]]==b,])
#   
#   # sample lower levels of hierarchy (if any)
#   if(length(cluster) > 1)
#     sub <- lapply(sub, .cases.resamp, cluster=cluster[-1], resample=resample[-1])
#   
#   # join and return samples
#   do.call(rbind, sub)
# }
# 
.cases.resamp <- function(model, .f, dat, cluster, resample) {
  # exit early for trivial data
  if(nrow(dat) == 1 || all(resample==FALSE))
    return(dat)
  
  # ver <- as.numeric_version(packageVersion("dplyr"))
  res <- dat
  
  for(i in 1:length(cluster)) {
    if(i==1 & resample[i]) {
      dots <- as.name(cluster[1])
      grouped <- dplyr::group_by(res, !!dots)
      g_rows <- dplyr::group_rows(grouped)
      # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- as.name(cluster[-i])
        grouped <- dplyr::group_by(res, .dots = !!dots) 
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- as.name(cluster[i])
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- purrr::map_dfr(res, function(df) { # ldply to purrr map from list to df
            grouped <- dplyr::group_by(df, .dots = !!dots)
            g_rows <- dplyr::group_rows(grouped)
            # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
            cls <- sample(seq_along(g_rows), replace = resample[i])
            idx <- unlist(g_rows[cls], recursive = FALSE)
            grouped[idx, ]
          }, .id = NULL)
        }
      }
    }
  }
  
  if(class(model) == "lmerMod"){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    tstar <- .f(lme4::lmer(formula = form, data = res, REML = reml)) 
    return(tstar)
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(class(model) == "lme"){
    tstar <- tryCatch(.f(updated.model(model = model, new.data = res)),  
                      error = function(e) e)
    if(inherits(tstar, "error")) {
      structure(rep(NA, length(.f(model))), fail.msgs = tstar$message)
    } else{
      tstar
    }
    return(tstar)
  } else if (class(model) == "glmerMod"){
    
  }else{
    stop("model class must be either 'lme' or 'lmerMod'")
  }
}