library(nlme)
library(boot)
library(RLRsim)

data(Socatt, package = "mlmRev")

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

# ==============================================================================
context("CGR bootstrap (lme)")
# ==============================================================================

mySumm <- function(.) { 
  c(beta = fixef(.), sigma = as.numeric(.$sigma), sig01 = as.numeric(VarCorr(.)[1,2]))
}

nsim <- 10

test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lme(mathAge11 ~ mathAge8 + gender + class,
                random = ~ 1 | school, data = jsp728)
  
  orig.stats <- mySumm(vcmodA)
  
  boo <- cgr_bootstrap(model = vcmodA, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "CGR")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lme(mathAge11 ~ mathAge8c + gender + class,
               random = ~ 1 | school, data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- cgr_bootstrap(model = rimod, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "CGR")
  expect_equal(boo$statistic, mySumm)
})


## See p. 34 of Goldstein's book
vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
              random = ~ 1 | school, data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- cgr_bootstrap(model = vcmodC, fn = mySumm, B = nsim)

test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "CGR")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lme(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class,
               random = ~ mathAge8c | school, data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- cgr_bootstrap(model = rcmod, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "CGR")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("three-level random intercept model",{
  skip_on_cran()
  rmA <- lme(rv ~ religion + year, random = ~ 1 | district/respond, data = Socatt)
  
  
  mySumm <- function(.) { 
    c(beta = fixef(.), sigma = as.numeric(.$sigma), sig01 = as.numeric(VarCorr(.)[1,2]))
  }
  
  orig.stats <- mySumm(rmA)
  
  
  boo <- cgr_bootstrap(model = rmA, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "CGR")
  expect_equal(boo$statistic, mySumm)
})