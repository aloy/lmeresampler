detach(package:lme4)
library(nlme)
library(nlmeU)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
set.seed(9221632)
trial <- parametric.lme(model = model, fn = fixef, B = 100)
trial

detach(package:nlme)
library(lme4)

# (model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
(model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))
set.seed(9221632)
trial2 <- parametric.lmerMod(model = model, fn = fixef, B = 10)
