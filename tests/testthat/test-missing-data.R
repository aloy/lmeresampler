skip_if_not_installed("glmmTMB")

library(lme4)
library(glmmTMB)

nsim <- 3

## incomplete cases: NA in the response (rows 10, 77, 120) and in the
## offset only (rows 4, 40)
set.seed(20261015)
dNA <- data.frame(f = factor(rep(1:30, each = 5)),
                  lo = log(runif(150, 1, 5)))
dNA$y <- rpois(150, exp(0.5 + dNA$lo + rnorm(30, sd = 0.5)[dNA$f]))
dNA$yg <- rnorm(150) + rnorm(30)[dNA$f]
dNA$y[c(10, 77, 120)] <- NA
dNA$yg[c(10, 77, 120)] <- NA
dNA$lo[c(4, 40)] <- NA

## glmmTMB's refit() looks up the data with eval.parent(), i.e. in the
## frame of its caller, so data local to this file are not found when
## refit() is called from inside lmeresampler; the data are placed in the
## global environment (and removed at the end of this file)
assign("dNA_lmeresampler_test", dNA, envir = globalenv())

## na.action is stored in the call and re-evaluated by refit(), so the
## models are fitted with literal na.action values rather than in a loop
mods <- list(
  lmer_omit = lmer(yg ~ 1 + (1 | f), data = dNA_lmeresampler_test),
  lmer_excl = lmer(yg ~ 1 + (1 | f), data = dNA_lmeresampler_test,
                   na.action = na.exclude),
  glmer_off = glmer(y ~ 1 + offset(lo) + (1 | f),
                    data = dNA_lmeresampler_test, family = poisson),
  tmb_omit = glmmTMB(yg ~ 1 + (1 | f), data = dNA_lmeresampler_test),
  tmb_excl = glmmTMB(yg ~ 1 + (1 | f), data = dNA_lmeresampler_test,
                     na.action = na.exclude),
  tmb_off = glmmTMB(y ~ 1 + offset(lo) + (1 | f),
                    data = dNA_lmeresampler_test, family = poisson)
)

## first fixed effect, for both lme4 and glmmTMB
fix1 <- function(m) unlist(fixef(m))[1]

test_that("refit_merMod aligns responses with incomplete cases", {
  for (nm in names(mods)) {
    m <- mods[[nm]]
    na.act <- attr(model.frame(m), "na.action")
    expect_length(na.act, if (grepl("off", nm)) 5 else 3)
    yobs <- model.response(model.frame(m))
    yfull <- rep(NA_real_, nrow(dNA))
    yfull[-na.act] <- yobs
    ## one value per observation used in the fit (residual bootstrap,
    ## glmmTMB's simulate()) or one per row of the data (lme4's simulate()
    ## with na.exclude): refitting the observed response must reproduce
    ## the original estimate in both cases
    for (y in list(yobs, yfull)) {
      res <- refit_merMod(data.frame(y), m, fix1)
      expect_equal(unname(res$tstar[[1]]), unname(fix1(m)),
                   tolerance = 1e-4, label = nm)
    }
  }
})

test_that("residual and parametric bootstraps work with incomplete cases", {
  skip_on_cran()
  for (nm in names(mods)) {
    m <- mods[[nm]]
    for (boot in list(resid_bootstrap, parametric_bootstrap)) {
      set.seed(20261015)
      boo <- boot(m, .f = extract_parameters, B = nsim)
      expect_s3_class(boo, "lmeresamp")
      expect_equal(nrow(boo$replicates), nsim, label = nm)
      expect_false(anyNA(boo$replicates), label = nm)
    }
  }
})

test_that("informative errors for bad .f output and failed refits", {
  m <- mods[["tmb_omit"]]
  expect_error(resid_bootstrap(m, .f = fixef, B = nsim),
               "function\\(m\\) fixef\\(m\\)\\$cond")
  ## .f succeeds on the original fit but fails on every refit
  f_refit_fails <- function(x) {
    if (identical(x, m)) fix1(x) else stop("boom")
  }
  expect_error(resid_bootstrap(m, .f = f_refit_fails, B = nsim),
               "all bootstrap refits failed; first error: boom")
})

rm("dNA_lmeresampler_test", envir = globalenv())
