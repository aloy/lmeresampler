library(lme4)
library(boot)
library(mlmRev)

##-------------------------------------------------------##
## Two-level models for JSP Data                         ##
##-------------------------------------------------------##


## Variance components model
## See p. 34 of Goldstein's book

vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + schoolMathAge8 + 
                 (1 | school), data = jsp728)

vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                 (1 | school), data = jsp728)

## Random coefficient model
## See p. 35 of Goldstein's book

rcmod <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class +  
                (mathAge8 | school), data = jsp728)


## Parametric bootstrap for the VC models

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

mySumm(vcmodA)

boo1 <- parametric_bootstrap.lmerMod(model = vcmodA, fn = mySumm, B = 10)

boo2 <- parametric_bootstrap.lmerMod(model = vcmodC, fn = mySumm, B = 10)


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
boo3 <- parametric_bootstrap.lmerMod(model = rmA, fn = mySumm2, B = 10)
