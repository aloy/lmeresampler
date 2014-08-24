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
t.res -> tstar


######
# WORKS UP TO HERE
######

tstar <- do.call("cbind", tstar) # Can these be nested?
rownames(tstar) <- names(t0)

RES <- structure(list(t0 = t0, t = t(tstar), R = B, data = model@frame,
                      seed = .Random.seed, statistic = fn,
                      sim = "parametric", call = match.call()),
                 class = "boot")

#### BAD OPT TEST #####
formula.mod <- formula(model)
formula.up <- update(old = formula.mod, new = ystar[,i] ~ .)
model.res <- lme(formula.up, data = sleepstudy, random = ~ re1|re2)


####### Singular test ######
model <- lme(ystar[,1] ~ Days, data = sleepstudy, random = ~Days|Subject)
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



