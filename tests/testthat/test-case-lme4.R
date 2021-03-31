library(lme4, quietly = TRUE)
library(dplyr, quietly = TRUE)

data(Socatt, package = "mlmRev")
Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal


# ==============================================================================
context("case bootstrap (lmerMod)")
# ==============================================================================
test_that("two-level additive random intercept model",{
  jjsp728 <- cbind(jsp728, .id = seq_len(nrow(jsp728)))
  grouped <- group_by(jjsp728, school) %>%
    summarise(count = n())
  
  cr1 <- .resamp.cases(dat = jjsp728, cluster = c("school", ".id"), resample = c(TRUE, TRUE))
  cr2 <- .resamp.cases(dat = jjsp728, cluster = c("school", ".id"), resample = c(FALSE, TRUE))
  cr3 <- .resamp.cases(dat = jjsp728, cluster = c("school", ".id"), resample = c(TRUE, FALSE))
  cr4 <- .resamp.cases(dat = jjsp728, cluster = c("school", ".id"), resample = c(FALSE, FALSE))
  
  expect_equal(nrow(cr2), nrow(jjsp728))
  expect_identical(cr4, jjsp728)
  expect_true(nrow(cr1) >= 48 * min(grouped$count))
  expect_true(nrow(cr1) <= 48 * max(grouped$count))
  expect_true(nrow(cr3) >= 48 * min(grouped$count))
  expect_true(nrow(cr3) <= 48 * max(grouped$count))
})

# ------------------------------------------------------------------------------

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

nsim <- 10

test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodA)
  
  boo <- case_bootstrap(model = vcmodA, .f = mySumm, B = nsim, resample = c(TRUE, TRUE))
  
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})




# ------------------------------------------------------------------------------
## See p. 97 of Goldstein's book
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- case_bootstrap(model = rimod, .f = mySumm, B = nsim, resample = c(TRUE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})

test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- case_bootstrap(model = vcmodC, .f = mySumm, B = nsim, resample = c(TRUE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- case_bootstrap(model = rcmod, .f = mySumm, B = nsim, resample = c(TRUE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("three-level random intercept model",{
  skip_on_cran()
  rmA <- lme4::lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
  
  orig.stats <- mySumm(rmA)
  boo <- case_bootstrap(model = rmA, .f = mySumm, B = nsim, resample = c(TRUE, TRUE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})


# ==============================================================================
context("case bootstrap (glmerMod)")
# ==============================================================================

mySumm <- function(.) { 
  c(beta = getME(., "beta"), sig01 = unname(getME(., "theta"))) 
}

test_that("two-level binomial logistic regression",{
  skip_on_cran()
  gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
  
  orig.stats <- mySumm(gm)
  boo <- case_bootstrap(model = gm, .f = mySumm, B = nsim, resample = c(TRUE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level poisson regression model",{
  skip_on_cran()
  gm <- glmer(TICKS ~ YEAR + cHEIGHT + (1|LOCATION),
              family="poisson", data=grouseticks)
  
  orig.stats <- mySumm(gm)
  boo <- case_bootstrap(model = gm, .f = mySumm, B = nsim, resample = c(TRUE, FALSE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})


test_that("three-level poisson regression model",{
  skip_on_cran()
  gm <- glmer(TICKS ~ YEAR + cHEIGHT + (1|LOCATION/BROOD),
              family="poisson",data=grouseticks)
  
  orig.stats <- mySumm(gm)
  boo <- case_bootstrap(model = gm, .f = mySumm, B = nsim, resample = c(FALSE, FALSE, TRUE))
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------