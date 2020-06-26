library(lme4)
library(lmeresampler)
library(bench)

## Package Examples

vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)
summary(vcmodA)

## you can write your own function to return stats, or use something like 'fixef'
mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

## running a parametric bootstrap 
set.seed(1234)
boo1 <- bootstrap(model = vcmodA, fn = mySumm, type = "parametric", B = 100, parallel = TRUE)

lb1 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "parametric", B = 10000, parallel = TRUE), filter_gc = FALSE)
lb1[c(2:9)]

## Not run: 
## running a cases bootstrap - only resampling the schools
boo2 <- bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE))

## running a cases bootstrap - resampling the schools and students within the school
boo2 <- bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE), parallel = FALSE)

lb2 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE), parallel = FALSE), filter_gc = FALSE)
lb2[c(2:9)]

## running a semi-parametric bootstrap
boo3 <- bootstrap(model = vcmodA, fn = mySumm, type = "cgr", B = 10000, parallel = TRUE)

lb3 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "cgr", B = 10000, parallel = TRUE), filter_gc = FALSE)
lb3[c(2:9)]

## running a residual bootstrap
boo4 <- bootstrap(model = vcmodA, fn = mySumm, type = "residual", B = 10000, parallel = TRUE)

lb4 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "residual", B = 10000, parallel = TRUE), filter_gc = FALSE)
lb4[c(2:9)]

## running an REB0 bootstrap
boo5 <- bootstrap(model = vcmodA, fn = mySumm, type = "reb", B = 100, reb_typ = 0)

lb5 <- bench::mark( bootstrap(model = vcmodA, fn = mySumm, type = "reb", B = 100, reb_typ = 0), filter_gc = FALSE)
lb5[c(2:9)]

