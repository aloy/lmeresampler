case.lmerMod <- function (model, fn, B, extra_step = FALSE){
  # TODO: put everything below into lapply to replicate
  .cases.resamp <- function (model, extra_step){
    # Draw sample of size J from level-2 units
    model.split <- split(x=model$data, f=model$data[3])
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
  
  # DEPRECATED rep.data <- as.data.frame( replicate(n = B, .cases.resamp(model = model, extra_step = extra_step)) )
  rep.data <- lapply(integer(B), eval.parent(substitute(function(...) .cases.resamp(model = model, extra_step = extra_step))))
  
  
  .cases.completion <- function(model, data, B, fn){
    t0 <- fn(model)
    
    # Refit the model and apply 'fn' to it using lapply
    form <- model@call$formula
    reml <- isREML(model)
    t.res <- matrix(0, ncol = 2, nrow = B)
    for(i in 1:B){
      myin <- ystar[,i]
      model.update <- update(object = model, fixed = myin ~ .)
      t.res[i,] <- fn(model.update)
    }
    t.res <- t(t.res)
    tstar <- split(t.res, rep(1:ncol(t.res), each = nrow(t.res)))
    
    tstar <- do.call("cbind", tstar) # Can these be nested?
    rownames(tstar) <- names(t0)
    
    RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                          seed = .Random.seed, statistic = fn,
                          sim = "parametric", call = match.call()),
                     class = "boot")
    
    return(RES)
  }
  
  # Plugin to .cases.completion due to small changes
  return(.cases.completion(model, rep.data, B, fn))
}
