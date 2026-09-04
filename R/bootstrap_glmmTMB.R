#' @rdname bootstrap
#' @export
#' @method bootstrap glmmTMB
bootstrap.glmmTMB <- function(model, .f = extract_parameters, type, B, resample,
                              reb_type, hccme,
                              aux.dist, orig_data = NULL, .refit = TRUE, rbootnoise = 0){

  if(type != "residual" && rbootnoise != 0) {
    stop("'rbootnoise' applicable only with residual bootstrapping. Do not define or use default 0.")
  }

  switch(type,
         parametric = parametric_bootstrap.glmmTMB(model, .f, B, .refit),
         residual = resid_bootstrap.glmmTMB(model, .f, B, .refit, rbootnoise),
         case = case_bootstrap.glmmTMB(model, .f, B, resample, orig_data, .refit),
         reb = stop("the REB bootstrap is not yet implemented for 'glmmTMB' objects."),
         wild = stop("the wild bootstrap is not available for 'glmmTMB' objects because hat values are not available."))
}


#' @rdname parametric_bootstrap
#' @export
#' @method parametric_bootstrap glmmTMB
#' @details
#' Identical to \code{\link{parametric_bootstrap.merMod}} -- neither
#' \code{simulate()} nor \code{\link{refit_merMod}} contain anything
#' lme4-specific, so the same implementation is reused as-is.
parametric_bootstrap.glmmTMB <- function(...) parametric_bootstrap.merMod(...)


#' @rdname case_bootstrap
#' @export
#' @method case_bootstrap glmmTMB
#' @details
#' Identical to \code{\link{case_bootstrap.merMod}}, which delegates the
#' class-specific pieces to the \code{.flist()} generic and to
#' \code{.resample_refit.cases()} (which uses \code{update()}, itself
#' generic over \code{merMod} and \code{glmmTMB}).
case_bootstrap.glmmTMB <- function(...) case_bootstrap.merMod(...)


#' @rdname resid_bootstrap
#' @export
#' @method resid_bootstrap glmmTMB
#' @details
#' Known limitation: for GLMMs, \code{weights(model)} returns \code{NULL}
#' for a \code{cbind(successes, failures)}-style binomial response (unlike
#' \code{glmer}, where it returns the trial counts), which makes the
#' response-resimulation step below error. Plain-Poisson and Gaussian
#' \code{glmmTMB} fits are unaffected.
resid_bootstrap.glmmTMB <- function(model, .f, B, .refit = TRUE, rbootnoise = 0){

  if(.refit) .f <- match.fun(.f)

  #Check the validity of rbootnoise
  if(!(rbootnoise >= 0 && rbootnoise <= 1)) {
    stop("'rbootnoise' between 0 to 1 should be used, such as 0.0001. The default 0 disables the feature of technical 2-level noise.")
  }

  setup <- .setup.glmmTMB(model, type = "residual", rbootnoise = rbootnoise)

  #For technical noise define the SD of e
  sde <- sd(setup[["e"]])

  #Calculate the number of clusters
  nclusters <- length(setup[["b"]][[1]][["(Intercept)"]])

  ystar <- as.data.frame(
    replicate(
      n = B,
      .resample.cgr(
        glmm = lme4::isGLMM(model),
        b = setup$b,
        e = setup$e,
        level.num = setup$level.num,
        Ztlist = setup$Ztlist,
        Xbeta = setup$Xbeta,
        vclist = setup$vclist,
        sig0 = setup$sig0,
        invlink = ifelse(lme4::isGLMM(model), stats::family(model)$linkinv, NULL),
        nclusters = nclusters,
        rbootnoise = rbootnoise,
        sde = sde
      )
    )
  )

  if(lme4::isGLMM(model)){
    fam <- stats::family(model)
    wts <- stats::weights(model)

    # simulate y
    simfun <- simfunList[[fam$family]]
    ystar <- purrr::map(ystar,
      ~simfun(model, nsim = 1, ftd = .x, wts = wts)
    )
  }

  if(!.refit) return(ystar)

  refits <- refit_merMod(ystar, model, .f)

  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "residual", warnings = refits$warnings)
}


#' @importFrom stats sigma
#' @export
#' @method extract_parameters glmmTMB
extract_parameters.glmmTMB <- function(model) {
  # VarCorr(model)$cond carries the same stddev/correlation/sc/useSc
  # attributes lme4 attaches, so reclassing lets us reuse lme4's own
  # as.data.frame method instead of duplicating it.
  vc <- glmmTMB::VarCorr(model)$cond
  class(vc) <- "VarCorr.merMod"
  vc <- as.data.frame(vc)

  c(
    beta = glmmTMB::fixef(model)$cond,
    vc = vc$vcov[is.na(vc$var2)]
  )
}


#' Determine whether a glmmTMB model is a GLMM (as opposed to an LMM)
#'
#' @description
#' A method for \code{lme4}'s \code{isGLMM()} generic. Registered onto
#' \code{lme4}'s method table in \code{.onLoad} (see \code{zzz.R}) rather
#' than via a static \code{S3method} directive, since \code{lme4} is only
#' an optional (\code{Suggests}) dependency of this package and importing
#' the generic would force it to load eagerly for every user.
#' @keywords internal
#' @noRd
isGLMM.glmmTMB <- function(x, ...) {
  fam <- stats::family(x)
  !(fam$family == "gaussian" && fam$link == "identity")
}
