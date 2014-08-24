library(nlme)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 10

## Begin Code

# Match function
fn <- match.fun(fn)
# Extract fixed effects
model.fixef <- fixed.effects(model)

# This works but do we want to import another package just to simulate?
library(nlmeU)
ystar <- simulateY(model, nsim = B)
row.names(ystar) <- 1:model$dims$N


## Where .bootstrap.completion starts
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


######
# WORKS UP TO HERE
######


RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")

