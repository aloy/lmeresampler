library(lme4)

data(Socatt, package = "mlmRev")
Socatt$religion <- relevel(Socatt$religion, ref = "none")
Socatt$rv <- as.numeric(as.character(Socatt$numpos))
Socatt$rv <- scale(Socatt$rv) # a plot shows this is clearly non-normal


# ==============================================================================
context("parametric bootstrap (lmerMod)")
# ==============================================================================

mySumm <- function(.) { 
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) 
}

nsim <- 10

jsp728$class <- relevel(jsp728$class, ref = "manual")


test_that("two-level additive random intercept model",{
  skip_on_cran()
  ## See p. 31 of Goldstein's book
  vcmodA <- lmer(mathAge11 ~ mathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodA)
  
  boo <- parametric_bootstrap.merMod(model = vcmodA, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------
test_that("two-level random intercept model without interaction",{
  skip_on_cran()
  ## See p. 97 of Goldstein's book
  rimod <- lmer(normAge11 ~ mathAge8c + gender + class + 
                  (1 | school), data = jsp728)
  
  orig.stats <- mySumm(rimod)
  boo <- parametric_bootstrap.merMod(model = rimod, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

test_that("two-level random intercept model with interaction",{
  skip_on_cran()
  ## See p. 34 of Goldstein's book
  vcmodC <- lmer(mathAge11 ~ mathAge8 * schoolMathAge8 + gender + class + 
                   (1 | school), data = jsp728)
  
  orig.stats <- mySumm(vcmodC)
  boo <- parametric_bootstrap.merMod(model = vcmodC, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("two-level random coefficient model with interaction",{
  skip_on_cran()
  ## See p. 35 of Goldstein's book
  rcmod <- lmer(mathAge11 ~ mathAge8c * schoolMathAge8 + gender + class + 
                  (mathAge8c | school), data = jsp728)
  
  orig.stats <- mySumm(rcmod)
  boo <- parametric_bootstrap.merMod(model = rcmod, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------

test_that("three-level random intercept model",{
  skip_on_cran()
  rmA <- lmer(rv ~ religion + year  + (1 | respond) + (1 | district), data = Socatt)
  
  orig.stats <- mySumm(rmA)
  boo <- parametric_bootstrap.merMod(model = rmA, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)   
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})


# model <- lme(mathgain ~ mathkind + sex + minority + ses, random = list( schoolid = ~mathkind, classid = ~1), 
# classroom, na.action = "na.omit")



# ==============================================================================
context("parametric bootstrap (glmerMod)")
# ==============================================================================

mySumm <- function(.) { 
  c(beta = getME(., "beta"), sig01 = unname(getME(., "theta"))) 
}

test_that("two-level binomial logistic regression",{
  skip_on_cran()
  gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
  
  orig.stats <- mySumm(gm)
  boo <- parametric_bootstrap(model = gm, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

# ------------------------------------------------------------------------------


test_that("two-level poisson regression model",{
  skip_on_cran()
  gm <- glmer(TICKS ~ YEAR + cHEIGHT + (1|LOCATION),
              family="poisson", data=grouseticks)
  
  orig.stats <- mySumm(gm)
  boo <- parametric_bootstrap(model = gm, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})


test_that("three-level poisson regression model",{
  skip_on_cran()
  gm <- glmer(TICKS ~ YEAR + cHEIGHT + (1|LOCATION/BROOD),
              family="poisson",data=grouseticks)
  
  orig.stats <- mySumm(gm)
  boo <- parametric_bootstrap(model = gm, .f = mySumm, B = nsim)
  
  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(unname(boo$stats$observed), unname(orig.stats))
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
  expect_equal(boo$.f, mySumm)
})

