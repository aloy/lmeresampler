#' Setup for resampling from lmerMod object
#'
#' @inheritParams bootstrap
#' @keywords internal
#' @noRd
.setup <- function(model, type, reb_type = NULL, rbootnoise = 0) {
  # Extract marginal means
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)

  # Extract ranef design matrix list
  Ztlist <- lme4::getME(object = model, name = "Ztlist")

  level.num <- lme4::getME(object = model, name = "n_rfacs")

  if (type == "reb") {
    Z <- lme4::getME(object = model, name = "Z")
    mresid <- lme4::getME(model, "y") - Xbeta

    # level-2 resid
    b <- solve(t(Z) %*% Z) %*% t(Z) %*% mresid # a single vector

    # level 1 resid
    e <- mresid - Z %*% b

    levs <- purrr::map(fl <- model@flist, levels)
    cnms <- lme4::getME(model, "cnms")
    b <- arrange_ranefs.lmerMod(b, fl, levs, cnms)
  }

  if (type == "residual") {
    # Extract and center random effects
    b <- purrr::map(lme4::ranef(model), .f = scale, scale = FALSE)
    b <- purrr::map(b, as.data.frame)

    # Extract and center error terms
    e <- scale(resid(model), scale = FALSE)
  }

  if (type == "wild") {
    mresid <- lme4::getME(model, "y") - Xbeta
    X <- lme4::getME(model, "X")
    .hatvalues <- diag(tcrossprod(X %*% solve(crossprod(X)), X))
    flist <- lme4::getME(model, "flist")
    if (length(flist) == 1) {
      n.lev <- nlevels(flist[[1]])
      flist <- as.numeric(flist[[1]])
    } else {
      flist <- lapply(flist, as.character)
      flist <- forcats::fct_inorder(do.call("paste", c(flist, sep = ":")))
      n.lev <- nlevels(flist)
      flist <- as.numeric(forcats::fct_inorder(flist))
    }
  }

  if (type == "residual" || type == "reb" && reb_type == 1) {
    sig0 <- stats::sigma(model)

    varcor <- lme4::VarCorr(model)

    #Add rbootnoise to 2-level variance when requested
    if (rbootnoise != 0) {
      varcor[[1]][[1]] <- varcor[[1]][[1]] +
        (attr(lme4::VarCorr(model), "sc") * rbootnoise)^2
    }
    vclist <- purrr::map(seq_along(b), .f = ~ bdiag(varcor[[names(b)[.x]]]))
    names(vclist) <- names(b)
  }

  if (type == "reb" && reb_type == 1) {
    if (level.num > 1) {
      stop("reb_type = 1 is not yet implemented for higher order models")
    }

    # Rescale u the residuals *prior* to resampling
    # so empirical variance is equal to estimated variance
    b <- purrr::map2(b, vclist, scale_center_ranef)
    e <- scale_center_e(e, sig0)
  }

  if (type %in% c("reb", "residual")) {
    RES <- list(Xbeta = Xbeta, b = b, e = e, Ztlist = Ztlist)

    if (type == "reb") {
      RES <- append(RES, list(flist = fl, levs = levs))
    } else {
      RES <- append(
        RES,
        list(level.num = level.num, sig0 = sig0, vclist = vclist)
      )
    }
  }

  if (type == "wild") {
    RES <- list(
      Xbeta = Xbeta,
      mresid = mresid,
      .hatvalues = .hatvalues,
      flist = flist,
      n.lev = n.lev
    )
  }

  RES
}


#' Setup for resampling from glmmTMB object
#'
#' @inheritParams bootstrap
#' @keywords internal
#' @noRd
.setup.glmmTMB <- function(model, type, rbootnoise = 0) {
  if (type != "residual") {
    stop("'.setup.glmmTMB' only supports the residual bootstrap.")
  }

  # Extract marginal means
  Xbeta <- predict(model, re.form = NA) # This is X %*% fixef(model)$cond

  reTrms <- .reTrms.glmmTMB(model)
  cnms <- reTrms$cnms
  Ztlist <- .make_Ztlist(cnms, reTrms$Gp, reTrms$Zt)

  level.num <- length(cnms)

  # Extract and center random effects
  b <- purrr::map(glmmTMB::ranef(model)$cond, .f = scale, scale = FALSE)
  b <- purrr::map(b, as.data.frame)

  # Extract and center error terms
  e <- scale(resid(model), scale = FALSE)

  sig0 <- stats::sigma(model)

  varcor <- glmmTMB::VarCorr(model)$cond

  #Add rbootnoise to 2-level variance when requested
  if (rbootnoise != 0) {
    varcor[[1]][[1]] <- varcor[[1]][[1]] + (attr(varcor, "sc") * rbootnoise)^2
  }
  vclist <- purrr::map(
    seq_along(b),
    .f = ~ Matrix::bdiag(varcor[[names(b)[.x]]])
  )
  names(vclist) <- names(b)

  list(
    Xbeta = Xbeta,
    b = b,
    e = e,
    Ztlist = Ztlist,
    level.num = level.num,
    sig0 = sig0,
    vclist = vclist
  )
}


#' Reconstruct the lme4-style random-effects structure for a glmmTMB model
#'
#' @description
#' glmmTMB's own \code{getME()} does not expose \code{Ztlist}/\code{Gp} the
#' way lme4's does (and its \code{"Gp"} entry currently errors due to an
#' upstream type bug). \code{cnms} and \code{flist}, however, are already
#' stored correctly -- computed once from the raw fitting data at fit time,
#' before any formula-side transformation -- in
#' \code{model$modelInfo$reTrms$cond}. \code{Gp} is then just the
#' cumulative per-term column count, and \code{Zt} is the transpose of
#' \code{getME(model, "Z")}, so the whole structure can be built without
#' ever re-evaluating the formula against data -- this matters because
#' \code{model.frame(model)} mangles some formula-side transformations into
#' a combined or renamed column (e.g. a \code{cbind(successes, failures)}
#' response, or \code{(log(x) | g)}), which re-evaluation would choke on.
#' @keywords internal
#' @noRd
.reTrms.glmmTMB <- function(model) {
  cnms <- model$modelInfo$reTrms$cond$cnms
  flist <- model$modelInfo$reTrms$cond$flist
  nc <- vapply(cnms, length, 1L)
  nl <- vapply(flist, nlevels, 1L)
  Gp <- c(0L, cumsum(nc * nl))
  Zt <- Matrix::t(glmmTMB::getME(model, "Z"))

  list(cnms = cnms, flist = flist, Gp = Gp, Zt = Zt)
}


#' Setup for resampling from lme object
#'
#' @inheritParams bootstrap
#' @keywords internal
#' @noRd
#' @importFrom HLMdiag extract_design
.setup.lme <- function(model, type, reb_type = NULL) {
  design <- HLMdiag::extract_design(model)

  # Extract marginal means
  Xbeta <- predict(model, level = 0)

  # Extract ranef design matrix list
  re.struct <- model$modelStruct$reStruct
  re.form <- formula(re.struct)

  Zlist <- extract_zlist.lme(model)

  level.num <- ncol(model$groups)

  if (type == "reb") {
    Z <- Matrix::Matrix(design$Z)
    mresid <- design$Y - Xbeta

    # level-2 resid
    b <- solve(t(Z) %*% Z) %*% t(Z) %*% mresid # a single vector

    # level 1 resid
    e <- mresid - Z %*% b

    levs <- purrr::map(fl <- model$groups, levels)
    cnms <- purrr::map(model$coefficients$random, colnames)
    b <- arrange_ranefs.lme(b, fl, levs, cnms)
  }
  if (type == "residual") {
    # Extract and center random effects
    reff <- nlme::ranef(model)
    if (level.num == 1) {
      reff <- list(reff)
      names(reff) <- colnames(model$groups)
    }

    b <- purrr::map(reff, .f = scale, scale = FALSE)
    b <- purrr::map(b, as.data.frame)

    # Extract and center error terms
    e <- scale(resid(model), scale = FALSE)
  }

  if (type == "wild") {
    mresid <- design$Y - Xbeta
    X <- design$X
    .hatvalues <- diag(tcrossprod(X %*% solve(crossprod(X)), X))
    flist <- model$groups
    flist <- flist[[ncol(flist)]]
    n.lev <- nlevels(flist)
    flist <- as.numeric(flist)
  }

  if (type == "residual" || type == "reb" && reb_type == 1) {
    sig0 <- stats::sigma(model)
    vclist <- purrr::map(
      as.matrix(model$modelStruct$reStruct),
      .f = ~ .x * sig0^2
    )
  }

  if (type == "reb" && reb_type == 1) {
    if (level.num > 1) {
      stop("reb_type = 1 is not yet implemented for higher order models")
    }

    # Rescale u the residuals *prior* to resampling
    # so empirical variance is equal to estimated variance
    b <- purrr::map2(b, vclist, scale_center_ranef)
    e <- scale_center_e(e, sig0)
  }

  if (type %in% c("reb", "residual")) {
    RES <- list(Xbeta = Xbeta, b = b, e = e, Zlist = Zlist)

    if (type == "reb") {
      RES <- append(RES, list(flist = fl, levs = levs))
    } else {
      RES <- append(
        RES,
        list(level.num = level.num, sig0 = sig0, vclist = vclist)
      )
    }
  }

  if (type == "wild") {
    RES <- list(
      Xbeta = Xbeta,
      mresid = mresid,
      .hatvalues = .hatvalues,
      flist = flist,
      n.lev = n.lev
    )
  }

  RES
}
