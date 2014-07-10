# Look at testdat package

# Data set
data(sleepstudy)
# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))


# Bootstrap function test
fm1.res <- bootstrap(model = fm1, fn = fixef, type = "par", B = 100)

# Parametric test
fm1.par.res <- parametric.lmerMod(model = fm1, fn = fixef, B = 100)
