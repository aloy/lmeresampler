library(lme4)
library(boot)
library(mlmRev)
library(dplyr)

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

# ==============================================================================
context("case bootstrap (lmerMod)")
# ==============================================================================

jjsp728 <- cbind(jsp728, .id = seq_len(nrow(jsp728)))
grouped <- group_by(jjsp728, school) %>%
  summarise(count = n())

cr1 <- .cases.resamp(dat = jjsp728, cluster = c("school", ".id"), replace = c(TRUE, TRUE))
cr2 <- .cases.resamp(dat = jjsp728, cluster = c("school", ".id"), replace = c(FALSE, TRUE))
cr3 <- .cases.resamp(dat = jjsp728, cluster = c("school", ".id"), replace = c(TRUE, FALSE))
cr4 <- .cases.resamp(dat = jjsp728, cluster = c("school", ".id"), replace = c(FALSE, FALSE))

cr1b <- .cases.resamp2(dat = jjsp728, cluster = c("school", ".id"), replace = c(TRUE, TRUE))
cr2b <- .cases.resamp2(dat = jjsp728, cluster = c("school", ".id"), replace = c(FALSE, TRUE))
cr3b <- .cases.resamp2(dat = jjsp728, cluster = c("school", ".id"), replace = c(TRUE, FALSE))
cr4b <- .cases.resamp2(dat = jjsp728, cluster = c("school", ".id"), replace = c(FALSE, FALSE))


test_that("two-level additive random intercept model",{
  expect_equal(nrow(cr2), nrow(jjsp728))
  expect_identical(cr4, jjsp728)
  expect_true(nrow(cr1) >= 48 * min(grouped$count))
  expect_true(nrow(cr1) <= 48 * max(grouped$count))
  expect_true(nrow(cr3) >= 48 * min(grouped$count))
  expect_true(nrow(cr3) <= 48 * max(grouped$count))
})

# ------------------------------------------------------------------------------

## See p. 31 of Goldstein's book
vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + 
                 (1 | school), data = jsp728)

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

orig.stats <- mySumm(vcmodA)

nsim <- 10

boo <- case_bootstrap(model = vcmodA, fn = mySumm, B = nsim, replace = c(TRUE, TRUE))

test_that("two-level additive random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})




# ------------------------------------------------------------------------------
## See p. 97 of Goldstein's book
rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                (1 | school), data = jsp728)

orig.stats <- mySumm(rimod)
boo <- case_bootstrap(model = rimod, fn = mySumm, B = nsim, replace = c(TRUE, TRUE))


test_that("two-level random intercept model without interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})


## See p. 34 of Goldstein's book
vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                 (1 | school), data = jsp728)

orig.stats <- mySumm(vcmodC)
boo <- case_bootstrap(model = vcmodC, fn = mySumm, B = nsim, replace = c(TRUE, TRUE))

test_that("two-level random intercept model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

## See p. 35 of Goldstein's book
rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                (mathAge8c | school), data = jsp728)

orig.stats <- mySumm(rcmod)
boo <- case_bootstrap(model = rcmod, fn = mySumm, B = nsim, replace = c(TRUE, TRUE))


test_that("two-level random coefficient model with interaction",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

rmA <- lme4::lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)

orig.stats <- mySumm(rmA)
boo <- case_bootstrap(model = rmA, fn = mySumm, B = nsim, replace = c(TRUE, TRUE, TRUE))


test_that("three-level random intercept model",{
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "parametric")
  expect_equal(boo$statistic, mySumm)
})