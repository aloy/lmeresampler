# Look at testdat package

# Data set
data(sleepstudy)

# fm1 model
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy))

# a user defined function
mySumm <- function(.) { 
  s <- sigma(.)
  c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}


# Bootstrap function test
fm1.res <- bootstrap(model = fm1, fn = fixef, type = "res", B = 100)
bootCI.mer(fm1.res, fm1, level = 0.95)

# Parametric test
fm1.par.res1 <- parametric.lmerMod(model = fm1, fn = fixef, B = 100)
fm1.par.res2 <- parametric.lmerMod(model = fm1, fn = mySumm, B = 100)

# Residual test
fm1.res.res1 <- residual.lmerMod(model = fm1, fn = fixef, B = 100)
fm1.res.res2 <- residual.lmerMod(model = fm1, fn = mySumm, B = 100)
