skip_if_not_installed("glmmTMB")

library(lme4)
library(glmmTMB)

data(sleepstudy, package = "lme4")
data(cbpp, package = "lme4")
data(grouseticks, package = "lme4")

nsim <- 3

# ==============================================================================
context("residual bootstrap (glmmTMB)")
# ==============================================================================

test_that("gaussian two-level random coefficient model", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  orig.stats <- extract_parameters(fm)
  set.seed(20260810)
  boo <- resid_bootstrap(model = fm, .f = extract_parameters, B = nsim)

  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(ncol(boo$replicates), length(orig.stats))
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "residual")
})

test_that("poisson two-level GLMM, via the bootstrap() dispatcher", {
  skip_on_cran()
  gm <- glmmTMB(TICKS ~ YEAR + cHEIGHT + (1 | LOCATION),
                family = "poisson", data = grouseticks)

  orig.stats <- extract_parameters(gm)
  set.seed(20260810)
  boo <- bootstrap(model = gm, .f = extract_parameters, type = "residual", B = nsim)

  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "residual")
})

# ==============================================================================
context("parametric bootstrap (glmmTMB)")
# ==============================================================================

test_that("gaussian two-level random coefficient model", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  orig.stats <- extract_parameters(fm)
  set.seed(20260810)
  boo <- parametric_bootstrap(model = fm, .f = extract_parameters, B = nsim)

  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "parametric")
})

# ==============================================================================
context("case bootstrap (glmmTMB)")
# ==============================================================================

test_that("gaussian two-level random coefficient model", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  orig.stats <- extract_parameters(fm)
  set.seed(20260810)
  boo <- case_bootstrap(model = fm, .f = extract_parameters, B = nsim,
                         resample = c(TRUE, TRUE))

  expect_equal(class(boo), "lmeresamp")
  expect_equal(boo$observed, orig.stats)
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
})

test_that("cbind() response requires orig_data (model.frame() mangles it otherwise)", {
  skip_on_cran()
  gm <- glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd),
                data = cbpp, family = binomial)

  set.seed(20260810)
  expect_warning(
    expect_error(
      case_bootstrap(model = gm, .f = extract_parameters, B = nsim,
                      resample = c(TRUE, TRUE)),
      "compatible sizes"
    ),
    "unnamed vectors"
  )

  set.seed(20260810)
  boo <- case_bootstrap(model = gm, .f = extract_parameters, B = nsim,
                         resample = c(TRUE, TRUE), orig_data = cbpp)

  expect_equal(class(boo), "lmeresamp")
  expect_equal(nrow(boo$replicates), nsim)
  expect_equal(boo$B, nsim)
  expect_equal(boo$type, "case")
})

# ==============================================================================
context("unsupported/unimplemented bootstraps (glmmTMB)")
# ==============================================================================

test_that("wild bootstrap is refused (no hat values available)", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  expect_error(
    bootstrap(model = fm, .f = extract_parameters, type = "wild", B = nsim,
              hccme = "hc2", aux.dist = "mammen"),
    "not available"
  )
})

test_that("reb bootstrap is not yet implemented", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  expect_error(
    bootstrap(model = fm, .f = extract_parameters, type = "reb", B = nsim, reb_type = 0),
    "not yet implemented"
  )
})

# ==============================================================================
context("glmmTMB utility methods")
# ==============================================================================

test_that("isGLMM.glmmTMB distinguishes LMMs from GLMMs", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  gm <- glmmTMB(TICKS ~ YEAR + cHEIGHT + (1 | LOCATION),
                family = "poisson", data = grouseticks)

  expect_false(lme4::isGLMM(fm))
  expect_true(lme4::isGLMM(gm))
})

test_that("extract_parameters.glmmTMB returns named beta/variance-component vector", {
  skip_on_cran()
  fm <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  params <- extract_parameters(fm)

  expect_type(params, "double")
  expect_true(all(c("beta.(Intercept)", "beta.Days") %in% names(params)))
  expect_equal(unname(params["beta.(Intercept)"]), unname(glmmTMB::fixef(fm)$cond["(Intercept)"]))
})
