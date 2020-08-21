## ----init, include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lmeresampler)

## ----fit models, results = FALSE, message = FALSE-----------------------------
library(lme4)
vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + (1 | school), data = jsp728)

library(nlme)
vcmodB <- lme(mathAge11 ~ mathAge8 + gender + class, random = ~1|school, data = jsp728)

## ----parameteric example------------------------------------------------------
# let's set .f = fixef to specify that we want only the fixed effects bootstrapped

# lme4
lmer_par_boot <- bootstrap(vcmodA, .f = fixef, type = "parametric", B = 100)

# nlme
lme_par_boot  <- bootstrap(vcmodB, .f = fixef, type = "parametric", B = 100)

## ----bootstrap output---------------------------------------------------------
# everything the lmeresamp object (and each of its sub-objects) contains
str(lmer_par_boot)

# an overview of those 10 list items
head(lmer_par_boot)

## ----residual example---------------------------------------------------------
# nlme
lme_res_boot <- bootstrap(vcmodB, .f = fixef, type = "residual", B = 100)

## ----cases example------------------------------------------------------------
# lme4
# resample the schools, but not the students
lmer_case_boot <- bootstrap(vcmodA, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))

# nlme
# do not resample the schools, but resample the students
lme_cases_boot1 <- bootstrap(vcmodB, .f = fixef, type = "case", B = 100, resample = c(FALSE, TRUE))

# nlme
# resample both the schools and the students
lme_cases_boot2 <- bootstrap(vcmodB, .f = fixef, type = "case", B = 100, resample = c(TRUE, TRUE))

## ----CGR example--------------------------------------------------------------
# lme4
lmer_cgr_boot <- bootstrap(vcmodA, .f = fixef, type = "cgr", B = 100)

# nlme
lme_cgr_boot  <- bootstrap(vcmodB, .f = fixef, type = "cgr", B = 100)

## ----REB example, message = FALSE, warning = FALSE----------------------------
# lme4
lmer_reb_boot0 <- bootstrap(vcmodA, .f = fixef, type = "reb", B = 100, reb_type = 0)
lmer_reb_boot1 <- bootstrap(vcmodA, .f = fixef, type = "reb", B = 100, reb_type = 1)

# nlme
lme_reb_boot2  <- bootstrap(vcmodB, .f = fixef, type = "reb", B = 100, reb_type = 2)

## ----print method, warning = FALSE, message = FALSE---------------------------
print(lmer_reb_boot0)

print(lme_reb_boot2, ci = TRUE) 

## ----confints, message = FALSE, warning = FALSE-------------------------------
confint(lme_reb_boot2) # all ci types, 0.95 confidence level

confint(lmer_reb_boot0, method = "boot-t", level= 0.90)

## ----parallel example, message = FALSE, warning = FALSE, cache = TRUE---------
library(purrr)
library(foreach)
library(doParallel)
set.seed(1234)
numCores <- 2

cl <- snow::makeSOCKcluster(numCores) # make a socket cluster
doParallel::registerDoParallel(cl) # how the CPU knows to run in parallel

b_parallel2 <- foreach(B = rep(250, 2), 
                       .combine = combine_lmeresamp,
                       .packages = c("lmeresampler", "lme4")) %dopar% {
  bootstrap(vcmodA, .f = fixef, type = "parametric", B = B)
}

stopCluster(cl)

## ----timings, message = FALSE, warning = FALSE, cache = TRUE------------------
library(tictoc)

tic()
b_nopar  <- bootstrap(vcmodA, .f = fixef, type = "parametric", B = 1000)
toc()

numCores <- 2
cl <- snow::makeSOCKcluster(numCores)
doParallel::registerDoParallel(cl) 

tic()
b_parallel2 <- foreach(B = rep(500, 2), 
                       .combine = combine_lmeresamp,
                       .packages = c("lmeresampler", "lme4", "purrr", "dplyr")) %dopar% {
  bootstrap(vcmodA, .f = fixef, type = "parametric", B = B)
}
toc()

stopCluster(cl) 

