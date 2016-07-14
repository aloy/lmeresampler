library(lme4)
library(boot)
library(mlmRev)

Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal

# ==============================================================================
context("REB bootstrap type = 0 (lmerMod)")
# ==============================================================================

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

nsim <- 10

test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lme4::lmer(mathAge11 ~ mathAge8 + gender + class + 
                         (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodA)
  
  boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- reb_bootstrap(model = rimod, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- reb_bootstrap(model = rcmod, fn = mySumm, B = nsim)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- try(reb_bootstrap(model = rmA, fn = mySumm, B = nsim))
# 
# 
# test_that("three-level random intercept model",{
#   expect_equal(class(boo), "try-error")
# })


# ==============================================================================
context("REB bootstrap type = 1 (lmerMod)")
# ==============================================================================

test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  mySumm <- function(.) { 
    s <- getME(., "sigma")
    c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
  }
  
  orig.stats <- mySumm(vcmodA)
  
  nsim <- 10
  
  boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- reb_bootstrap(model = rimod, fn = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})


test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- reb_bootstrap(model = rcmod, fn = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- reb_bootstrap(model = rmA, fn = mySumm, B = nsim, reb_type = 1)
# 
# 
# test_that("three-level random intercept model",{
#   expect_equal(class(boo), "boot")
#   expect_equal(boo$t0, orig.stats)
#   expect_equal(nrow(boo$t), nsim)
#   expect_equal(ncol(boo$t), length(orig.stats))
#   expect_equal(boo$R, nsim)
#   expect_equal(boo$sim, "reb")
#   expect_equal(boo$statistic, mySumm)
# })

# ==============================================================================
context("REB bootstrap type = 2 (lmerMod)")
# ==============================================================================

mySumm <- function(.) {
  c(beta = lme4::fixef(.), sigma =c(diag(bdiag(lme4::VarCorr(.))), lme4::getME(., "sigma")^2))
}

test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodA)
  
  nsim <- 10
  
  boo <- reb_bootstrap(model = vcmodA, fn = mySumm, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- reb_bootstrap(model = rimod, fn = mySumm, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})


test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, fn = mySumm, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- reb_bootstrap(model = rcmod, fn = mySumm, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "boot")
  expect_equal(boo$t0, orig.stats)
  expect_equal(nrow(boo$t), nsim)
  expect_equal(ncol(boo$t), length(orig.stats))
  expect_equal(boo$R, nsim)
  expect_equal(boo$sim, "reb")
  expect_equal(boo$statistic, mySumm)
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- reb_bootstrap(model = rmA, fn = mySumm, B = nsim, reb_type = 2)
# 
# 
# test_that("three-level random intercept model",{
#   expect_equal(class(boo), "boot")
#   expect_equal(boo$t0, orig.stats)
#   expect_equal(nrow(boo$t), nsim)
#   expect_equal(ncol(boo$t), length(orig.stats))
#   expect_equal(boo$R, nsim)
#   expect_equal(boo$sim, "reb")
#   expect_equal(boo$statistic, mySumm)
# })