parametric.lme <- function(model, fn, B){
  # Match function
  fn <- match.fun(fn)
  # Extract fixed effects
  model.fixef <- fixed.effects(model)
  
  ystar <- simulateY(model, nsim = B)
  row.names(ystar) <- 1:model$dims$N
  
  t0 <- fn(model)
  
  t.res <- list()
  length(t.res) <- B
#    t.res <- matrix(0, ncol = 2, nrow = B)
#   for(i in 1:B){
#     myin <- ystar[,i]
#     model.update <- update(object = model, fixed = myin ~ .)
#     t.res[i,] <- fn(model.update)
#   }

  for(i in 1:B){
#     myin <- ystar[,i]
#     model.update <- nlme:::update.lme(object = model, fixed = myin ~ .)
#     t.res[i,] <- fn(model.update)
    try.fit <- try( updated.model(model = model, new.y = ystar[,i]) )
    ### NOTE: Need to check if there was an error in try, and then use fn()
    t.res[[i]] <- updated.model(model = model, new.y = ystar[,i])
    
  }
  t.res <- do.call('rbind', t.res)
  tstar <- data.frame(t(t.res))

#   tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
#   
#   tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  colnames(tstar) <- paste("sim", 1:ncol(tstar), sep = "_")
  
  RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model$data,
                        seed = .Random.seed, statistic = fn,
                        sim = "parametric", call = match.call()),
                   class = "boot")
  return(RES)
}
