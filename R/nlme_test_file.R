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


t.res <- matrix(0, ncol = ncol(ystar), nrow = nrow(ystar))
for(i in 1:B){ystar[,i]
  myin <- ystar[,i]
  model.update <- lme.refit(model, myin)
  t.res[,i] <- fn(model.update)
}

lme.refit <- function(model, fixed.update){
  res <- update(object = model, fixed = fixed.update ~ .)
  return(res)
}

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")




####### Singular test ######
model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 1

### BEGIN PAR CODE ###
fn <- match.fun(fn)

model.fixef <- fixed.effects(model) # Extract fixed effects
ystar <- simulate.lme.data(model, nsim = B, na.action = na.exclude)

t0 <- fn(model)

# Refit the model and apply 'fn' to it using lapply
model.update <- update(object = model, ystar ~ .)
anova(model.update, model) # This should not run if the update worked
tstar <- fn(model.update)

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")



