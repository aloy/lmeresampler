# Look at testdat package

# Data set
data(sleepstudy)
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))
parametric.lmerMod(model = fm1, fn = fixef, B = 1000)