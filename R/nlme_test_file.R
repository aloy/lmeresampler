library(nlme)
library(nlmeU)

data(sleepstudy, package = "lme4")

model <- lme(Reaction ~ Days, data = sleepstudy, random = ~Days|Subject)
fn <- fixed.effects
B <- 10


library(lme4)
(model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
#(model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))
B <- 10
fn <- fixef