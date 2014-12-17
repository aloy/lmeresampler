# Look at testdat package

library(lme4)
library(boot)
library(HLMdiag)

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
#bootCI.mer(fm1.res, fm1, level = 0.95)

# Parametric test
fm1.par.res1 <- parametric.lmerMod(model = fm1, fn = fixef, B = 100)
fm1.par.res2 <- parametric.lmerMod(model = fm1, fn = mySumm, B = 100)

# Residual test
fm1.res.res1 <- residual.lmerMod(model = fm1, fn = fixef, B = 100)
fm1.res.res2 <- residual.lmerMod(model = fm1, fn = mySumm, B = 100)

# Cases test
fm1.case.res1 <- case.lmerMod(model = fm1, fn = fixef, B = 100)
fm1.case.res2 <- case.lmerMod(model = fm1, fn = mySumm, B = 100)
fm1.case.e.res1 <- case.lmerMod(model = fm1, fn = fixef, B = 100, extra_step = TRUE)
fm1.case.e.res2 <- case.lmerMod(model = fm1, fn = mySumm, B = 100, extra_step = TRUE)

# CGR
fm1.cgr.res1 <- cgr.lmerMod(model = fm1, fn = fixef, B = 100)



# Looking at a three level model
library(WWGbook)
data(classroom)

fm2 <- lmer(mathgain ~ mathkind + sex + minority + ses + housepov + (1|schoolid) + (1|classid), 
            data = classroom, na.action = "na.omit")

fm3 <- lmer(mathgain ~ mathkind + sex + minority + ses + housepov + (ses|schoolid) + (1|classid), 
            data = classroom, na.action = "na.omit")


# Parametric test
fm2.par.res1 <- parametric.lmerMod(model = fm2, fn = fixef, B = 100)
fm2.par.res2 <- parametric.lmerMod(model = fm2, fn = mySumm, B = 100)

# Residual test
fm2.res.res1 <- residual.lmerMod(model = fm2, fn = fixef, B = 100)
fm2.res.res2 <- residual.lmerMod(model = fm2, fn = mySumm, B = 100)

fm3.res.res1 <- residual.lmerMod(model = fm3, fn = fixef, B = 100)
fm3.res.res2 <- residual.lmerMod(model = fm3, fn = mySumm, B = 100)

# Cases test
fm2.case.res1 <- case.lmerMod(model = fm2, fn = fixef, B = 100)
fm2.case.res2 <- case.lmerMod(model = fm2, fn = mySumm, B = 100)
fm2.case.e.res1 <- case.lmerMod(model = fm2, fn = fixef, B = 100, extra_step = TRUE)
fm2.case.e.res2 <- case.lmerMod(model = fm2, fn = mySumm, B = 100, extra_step = TRUE)


# Looking at model from Goldstein
library(mlmRev)
Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

mod <- lmer(rv ~ religion + year + (1 | respond) + (1 | district), data = Socatt)
bse <- parametric.lmerMod(model = mod, fn = fixef, B = 1000)
