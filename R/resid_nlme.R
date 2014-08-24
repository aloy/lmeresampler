library(nlme)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 10

## Begin Code

fn <- match.fun(fn)

ystar <- as.data.frame( replicate(n = B, .resample.resids(model = model)) )

.resample.resids <- function(model){
  
  # Extract fixed part of the model
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)
  
  # Extract random effects
  model.ranef <- random.effects(model)
  
  # Extract residuals
  model.resid <- residuals(model)
  
  # Extract Z design matrix
  # Z <- getME(object = model, name = "Ztlist") # FIX THIS

  J <- nrow(model.ranef)
  
  # Create Z design matrix
  Z <- diag(J)
                    
  # Sample of b*
  bstar.index <- sample(x = seq_len(J), size = J, replace = TRUE)
  bstar <- model.ranef[bstar.index,]
  
  level.num <- model$dims$Q # I believe this is where it lists it?
  
  if(level.num == 1){
    bstar <- lapply(bstar, FUN = function(x) as.list(x))[[1]]
    #names(bstar) <- names(Z)
  } else {
    bstar <- lapply(bstar, FUN = function(x) as.data.frame(x))
    bstar <- do.call(c, bstar)
    #names(bstar) <- names(Z)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = bstar, zstar = Z)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  
  # Resample residuals
  estar <- sample(x = model.resid, size = length(model.resid), replace = TRUE)
  
  # Combine function
  y.star <- as.numeric(Xbeta + Zbstar.sum + estar)
  
  return(y.star)
  
}

t0 <- fn(model)

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