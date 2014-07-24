extra_step <- FALSE

bootstrap.completion <- function(model, ystar, B, fn){
  t0 <- fn(model)
  
  # Refit the model and apply 'fn' to it using lapply
  tstar <- lapply(ystar, function(x) {
    fn(refit(object = model, newresp = x))
  })
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  
  return(RES)
}

cases.completion <- function(model, data, B, fn){
  t0 <- fn(model)
  
  # Refit the model and apply 'fn' to it using lapply
  tstar <- lapply(data, function(x) {
    fn(lmer(lmer(model@call$formula, data, isREML(model))))
  })
  
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  
  return(RES)
}

data(sleepstudy)

# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))

cases.resamp <- function (model, extra_step){
  model.split <- split(x=model@frame, f=model@flist)
  model.split.samp <- sample(x=model.split, size = length(model.split), replace = TRUE)
  # For each sample, draw a sample of the cases from the level-2 unit
  
  if(extra_step == TRUE){
    model.resamp <- lapply(model.split.samp,
                         FUN = function(x) {
                           J <- nrow(x)
                           
                           # Sample of level-2 rows
                           model.sub.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                           resampled <- x[model.sub.index,]
                           return(resampled)
                         })
    model.comb <- do.call('rbind', model.resamp)
  } else{ # else statement needs to be located here
    model.comb <- do.call('rbind', model.split.samp)
  }
}

B <- 100
ystar <- as.data.frame( replicate(n = B, cases.resamp(model = fm1, extra_step = extra_step)) )

#' ISS: I think I am passing in ystar as a list of simulations but then
#' it throws an error:
#' Error in refit.merMod(object = model, newresp = x) : 
#' refit not implemented for lists with length>1: consider ‘lapply(object,refit)’
#' which does not make sense because that is what we do with all the other
#' bootstraps above
fm1.RES <- bootstrap.completion(fm1, ystar, B = B, fn = fixef)

#####
#' Noticable issues:
#' this outputs a list in the proper format but
#' it does not give the proper results. I have tried to follow steps and figure
#' out what the issue is but have not found anything at this point. I believe it
#' is only happening for one simulation