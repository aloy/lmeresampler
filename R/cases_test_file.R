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

data(sleepstudy)

# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))

fm1.split <- split(x=fm1@frame, f=fm1@flist)
fm1.split.samp <- sample(x=fm1.split, size = length(fm1.split), replace = TRUE)
# For each sample, draw a sample of the cases from the level-2 unit

if(extra_step == TRUE){
  fm1.resamp <- lapply(fm1.split.samp,
                         FUN = function(x) {
                           J <- nrow(x)
                           
                           # Sample of level-2 rows
                           fm1.sub.index <- sample(x = seq_len(J), size = J, replace = TRUE)
                           resampled <- x[fm1.sub.index,]
                           return(resampled)
                         })
  fm1.comb <- do.call('rbind', fm1.resamp)
}
else{
  fm1.comb <- do.call('rbind', fm1.split.samp)
}

fm1.RES <- bootstrap.completion(fm1, fm1.comb, B = 100, fn = fixef)