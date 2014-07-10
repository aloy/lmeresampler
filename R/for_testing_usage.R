# Look at testdat package

# Data set
data(sleepstudy)
# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))

# Parametric test
parametric.lmerMod(model = fm1, fn = fixef, B = 100)
