library(nlme)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 10


# This works but do we want to import another package just to simulate?
library(nlmeU)
sims <- simulateY(model, nsim = B)

### BEGIN PAR CODE ###
fn <- match.fun(fn)

model.fixef <- fixed.effects(model) # Extract fixed effects


t0 <- fn(model)


# Originally I had this: 
# tstar <- lapply(ystar, function(x) {
#   update(object = model, x ~ .)
#   fn(model)
# })
# But ran into this error:
# Error in eval(expr, envir, enclos) : object 'x' not found
# And when I change it to Reaction it runs but I do not think it works
# Refit the model and apply 'fn' to it using lapply
tstar <- lapply(ystar, function(x) {
  model.update <- update(object = model, fixed = x ~ .)
  t.res <- fn(model.update)
  return(t.res)
})

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



