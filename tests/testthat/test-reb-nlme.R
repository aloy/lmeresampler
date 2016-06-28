library(nlme, quietly = TRUE)
library(boot)
library(RLRsim)

data(Socatt, package = "mlmRev")

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

# ==============================================================================
context("REB bootstrap type = 0 (lme)")
# ==============================================================================

## See p. 31 of Goldstein's book
vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class,
              random = ~ 1 | school, data = jsp728)


mySumm <- function(.) { 
  c(beta = fixef(.), sigma = as.numeric(.$sigma), sig01 = as.numeric(VarCorr(.)[1,2]))
}

orig.stats <- mySumm(vcmodA)

nsim <- 10

set.seed(7142015)
boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim, reb_type = 0)

test_that("two-level additive random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb0")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 34 of Goldstein's book
vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
              random = ~ 1 | school, data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim)

test_that("two-level random intercept model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb0")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 35 of Goldstein's book
rcmod <- lme(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class,
             random = ~ mathAge8c | school, data = jsp728)

orig.stats <- mySumm(rcmod)
boo <- reb_bootstrap(model = rcmod, fn = mySumm, B = nsim)


test_that("two-level random coefficient model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb0")
  expect_equal(boo$statistic, mySumm)
})


# ==============================================================================
context("REB bootstrap type = 1 (lme)")
# ==============================================================================


## See p. 31 of Goldstein's book
vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class,
              random = ~ 1 | school, data = jsp728)


mySumm <- function(.) { 
  c(beta = fixef(.), sigma = as.numeric(.$sigma), sig01 = as.numeric(VarCorr(.)[1,2]))
}

orig.stats <- mySumm(vcmodA)

nsim <- 10

set.seed(7142015)
boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim, reb_type = 1)

test_that("two-level additive random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb1")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 34 of Goldstein's book
vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
              random = ~ 1 | school, data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim, reb_type = 1)

test_that("two-level random intercept model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb1")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 35 of Goldstein's book
rcmod <- lme(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class,
             random = ~ mathAge8c | school, data = jsp728)

orig.stats <- mySumm(rcmod)
boo <- reb_bootstrap(model = rcmod, fn = mySumm, B = nsim, reb_type = 1)


test_that("two-level random coefficient model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb1")
  expect_equal(boo$statistic, mySumm)
})


# ==============================================================================
context("REB bootstrap type = 2 (lme)")
# ==============================================================================



## See p. 31 of Goldstein's book
vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class,
              random = ~ 1 | school, data = jsp728)


mySumm <- function(.) {
  c(beta = nlme::fixef(.), sigma = c(diag(nlme::getVarCov(.)), .$sigma^2))
}

orig.stats <- mySumm(vcmodA)

nsim <- 10

set.seed(7142015)
boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim, reb_type = 2)

test_that("two-level additive random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb2")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 34 of Goldstein's book
vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
              random = ~ 1 | school, data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim, reb_type = 2)

test_that("two-level random intercept model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb2")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 35 of Goldstein's book
rcmod <- lme(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class,
             random = ~ mathAge8c | school, data = jsp728)

orig.stats <- mySumm(rcmod)
boo <- reb_bootstrap(model = rcmod, B = nsim, reb_type = 2)


test_that("two-level random coefficient model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb2")
  expect_equal(boo$statistic, mySumm)
})
