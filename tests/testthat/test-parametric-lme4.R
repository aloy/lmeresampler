library(lme4)
library(boot)
library(mlmRev)

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal


context("parametric bootstrap (lmerMod)")


vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + schoolMathAge8 + 
                 (1 | school), data = jsp728)

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

orig.stats <- mySumm(vcmodA)

nsim <- 10

boo <- parametric_bootstrap.lmerMod(model = vcmodA, fn = mySumm, B = nsim)

test_that("two-level additive random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})


vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                 (1 | school), data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- parametric_bootstrap.lmerMod(model = vcmodC, fn = mySumm, B = nsim)

test_that("two-level random intercept model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})




rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)

orig.stats <- mySumm(rmA)
boo <- parametric_bootstrap.lmerMod(model = rmA, fn = mySumm, B = nsim)


test_that("three-level random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})