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
boo1 <- bootstrap(model = vcmodA, fn = mySumm, type = "parametric", B = 100, parallel = TRUE, nCores = 2)

lb1 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "parametric", B = 100, parallel = TRUE), filter_gc = FALSE)
lb1[c(2:9)]

## Not run: 
## running a cases bootstrap - only resampling the schools
boo2 <- bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE))

## running a cases bootstrap - resampling the schools and students within the school
boo2 <- bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE))

lb2 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "case", B = 100, resample = c(TRUE, FALSE)), filter_gc = FALSE)
lb2[c(2:9)]

## running a semi-parametric bootstrap
boo3 <- bootstrap(model = vcmodA, fn = mySumm, type = "cgr", B = 100)

lb3 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "cgr", B = 100), filter_gc = FALSE)
lb3[c(2:9)]

## running a residual bootstrap
boo4 <- bootstrap(model = vcmodA, fn = mySumm, type = "residual", B = 100)

lb4 <- bench::mark(bootstrap(model = vcmodA, fn = mySumm, type = "residual", B = 100), filter_gc = FALSE)
lb4[c(2:9)]

## running an REB0 bootstrap
boo5 <- bootstrap(model = vcmodA, fn = mySumm, type = "reb", B = 100, reb_typ = 0)

lb5 <- bench::mark( bootstrap(model = vcmodA, fn = mySumm, type = "reb", B = 100, reb_typ = 0), filter_gc = FALSE)
lb5[c(2:9)]

# nCores idea, this is in generics.R, easier for me to see where I made changes here

bootstrap <- function(model, fn, type, B, resample = NULL, reb_type = NULL, parallel = FALSE, nCores = NULL) {
  if(!type %in% c("parametric", "residual", "case", "cgr", "reb"))
    stop("'type' must be one of 'parametric', 'residual', 'case', 'cgr', or 'reb'")
  if(!is.null(reb_type))
    if(!reb_type %in% 0:2) 
      stop("'reb_type' must be either 0, 1, or 2")
  if(parallel == FALSE) nCores <- 1
  else {
    if(nCores %in% 2:parallel::detectCores())
      stop("for parallelization 'nCores' must be greater than 1 and within the range of your machine's cores")
    if(is.null(nCores)) nCores <- 2
  }
  UseMethod("bootstrap", model)
}
