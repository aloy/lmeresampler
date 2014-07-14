# Look at testdat package

# Data set
data(sleepstudy)

# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))


# Bootstrap function test
fm1.res <- bootstrap(model = fm1, fn = fixef, type = "res", B = 100)
bootCI.mer(fm1.res, fm1, level = 0.95)

# Parametric test
fm1.par.res <- parametric.lmerMod(model = fm1, fn = fixef, B = 100)

# Residual test
fm1.res.res <- residual.lmerMod(model = fm1, fn = fixef, B = 100)
