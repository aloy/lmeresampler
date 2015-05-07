library(nlme)
library(boot)
library(nlmeU) # for simulateY

##-------------------------------------------------------##
## Two-level models for JSP Data                         ##
##-------------------------------------------------------##


## Variance components model
## See p. 34 of Goldstein's book

vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class + schoolMathAge8, 
              data = jsp728, random = ~ 1 | school)

vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
              data = jsp728, random = ~ 1 | school)

## Random coefficient model
## See p. 35 of Goldstein's book

rcmod <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
             data = jsp728, random = ~ mathAge8 | school)


## Parametric bootstrap for the VC models

fixef(vcmodA)

boo1 <- parametric.lme(model = vcmodA, fn = fixef, B = 100)

boo2 <- parametric.lmerMod(model = vcmodC, fn = mySumm, B = 1000)


##-------------------------------------------------------##
## Three-level models for 
##-------------------------------------------------------##

## Social attitidues survey
## See p. 87 and 97 of Goldstein's book

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)

mySumm2 <- function(.) { 
  s <- sigma(.)
  c(fixef(.), sigma = s, sig = unname(s * getME(., "theta"))) 
}

mySumm2(rmA)
boo3 <- parametric.lmerMod(model = rmA, fn = mySumm2, B = 1000)
