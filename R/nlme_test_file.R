library(nlme)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixef
B <- 10

### BEGIN PAR CODE ###
fn <- match.fun(fn)

model.fixef <- fixed.effects(model) # Extract fixed effects
ystar <- simulate.lme.data(model, nsim = B, na.action = na.exclude)

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
tstar <- lapply(ystar, function(Reaction) {
  update(object = model, Reaction ~ .)
  fn(model)
})

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")