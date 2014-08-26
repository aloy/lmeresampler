parametric.lme <- function(model, fn, B){
  # Match function
  fn <- match.fun(fn)
  # Extract fixed effects
  model.fixef <- fixed.effects(model)
  
  ystar <- simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  
  t0 <- fn(model)
  
   t.res <- matrix(0, ncol = 2, nrow = B)
#   for(i in 1:B){
#     myin <- ystar[,i]
#     model.update <- update(object = model, fixed = myin ~ .)
#     t.res[i,] <- fn(model.update)
#   }

  for(i in 1:B){
    t.res[i,] <- updated.model(model = model, up.reaction = ystar[,i])
  }
  tstar <- data.frame(t(t.res))

#   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
#   
#   tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  colnames(tstar) <- 1:ncol(tstar)
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
}


# Currently getting an error again...
# Maybe use the 'reformulate' fn
updated.model<- function(model, up.reaction){
  # Extract formulas and data
  mod.fixd <- as.formula(model$call$fixed)
  mod.rand <- as.formula(model$call$random)
  mod.data <- model$data
  # Place ystars in data
  mod.data$Reaction <- up.reaction
  # create new lme
  out.lme <- lme(fixed = mod.fixd, data = mod.data, random = mod.rand)
  return(out.lme)
}