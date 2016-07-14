library(nlme, quietly = TRUE)
library(boot)
library(RLRsim)

data(Socatt, package = "mlmRev")

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

# ==============================================================================
context("residual bootstrap (lme)")
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
  
  set.seed(7142015)
  boo <- resid_bootstrap(model = vcmodA, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "resid")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lme(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class, 
                random = ~ 1 | school, data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- resid_bootstrap(model = vcmodC, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "resid")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lme(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class,
               random = ~ mathAge8c | school, data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- resid_bootstrap(model = rcmod, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "resid")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("three-level random intercept model",{
  skip_on_cran()
  rmA <- lme(rv ~ religion + year, random = ~ 1 | district/respond, data = Socatt)
  
  
  mySumm <- function(.) { 
    c(beta = fixef(.), sigma = as.numeric(.$sigma), 
      sig.dist = as.numeric(VarCorr(.)[2,2]), 
      sig.int = as.numeric(VarCorr(.)[4,2]))
  }
  
  orig.stats <- mySumm(rmA)
  
  boo <- resid_bootstrap(model = rmA, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "resid")
  expect_equal(boo$statistic, mySumm)
})

