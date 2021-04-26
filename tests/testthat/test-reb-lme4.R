library(lme4)

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
  
  boo <- reb_bootstrap(model = vcmodA, .f = mySumm, B = nsim, reb_type = 0)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb0")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- reb_bootstrap(model = rimod, .f = mySumm, B = nsim, reb_type = 0)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb0")
  expect_equal(boo$.f, mySumm)
})

test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, .f = mySumm, B = nsim, reb_type = 0)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb0")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- reb_bootstrap(model = rcmod, .f = mySumm, B = nsim, reb_type = 0)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb0")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- try(reb_bootstrap(model = rmA, .f = mySumm, B = nsim))
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
  
  boo <- reb_bootstrap(model = vcmodA, .f = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb1")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- reb_bootstrap(model = rimod, .f = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb1")
  expect_equal(boo$.f, mySumm)
})


test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, .f = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb1")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- reb_bootstrap(model = rcmod, .f = mySumm, B = nsim, reb_type = 1)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb1")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- reb_bootstrap(model = rmA, .f = mySumm, B = nsim, reb_type = 1)
# 
# 
# test_that("three-level random intercept model",{
#   expect_equal(class(boo), "boot")
#   expect_equal(boo$t0, orig.stats)
#   expect_equal(nrow(boo$t), nsim)
#   expect_equal(ncol(boo$t), length(orig.stats))
#   expect_equal(boo$B, nsim)
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
  
  orig.stats <- extract_parameters(vcmodA)
  
  nsim <- 10
  
  boo <- reb_bootstrap(model = vcmodA, .f = extract_parameters, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb2")
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- extract_parameters(rimod)
  boo <- reb_bootstrap(model = rimod, .f = extract_parameters, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb2")
})


test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- extract_parameters(vcmodC)
  boo <- reb_bootstrap(model = vcmodC, .f = mySumm, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb2")
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- extract_parameters(rcmod)
  boo <- reb_bootstrap(model = rcmod, .f = extract_parameters, B = nsim, reb_type = 2)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "reb2")
})

# ------------------------------------------------------------------------------

# rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
# 
# orig.stats <- mySumm(rmA)
# boo <- reb_bootstrap(model = rmA, .f = mySumm, B = nsim, reb_type = 2)
# 
# 
# test_that("three-level random intercept model",{
#   expect_equal(class(boo), "boot")
#   expect_equal(boo$t0, orig.stats)
#   expect_equal(nrow(boo$t), nsim)
#   expect_equal(ncol(boo$t), length(orig.stats))
#   expect_equal(boo$B, nsim)
#   expect_equal(boo$sim, "reb2")
#   expect_equal(boo$statistic, mySumm)
# })