scale_center_ranef <- function(b, vc){
  u <- scale(b, scale = FALSE)
  S <- (t(u) %*% u) / nrow(u)
  R <- vc
  
  Ls <- t(chol(S, pivot = TRUE))
  Lr <- t(chol(R, pivot = TRUE))
  A <- t(Lr %*% solve(Ls))
  
  Uhat <- as.matrix(u %*% A)
  data.frame(Uhat)
}

scale_center_e <- function(x, sigma) {
  estar <- sigma * x %*% ((t(x) %*% x) / length(x))^(-1/2)
  scale(estar, scale = FALSE)
}


extract_parameters.lmerMod <- function(model) {
  sig.e <- sigma(model)
  vc <- as.data.frame(lme4::VarCorr(model))
  
  c(
    beta = getME(model, "beta"), 
    vc = vc$vcov[is.na(vc$var2)]
  )
}


#' Organize ranef vector into lists
#' 
#' @param b ranef estimates
#' @param fl flist from lmerMod object, a list of the grouping 
#'   variables (factors) involved in the random effect terms
#' @param levs unique levels of each factor in fl
#' @param cnms a list component names for each ranef
#' @keywords internal
#' @noRd
arrange_ranefs.lmerMod <- function(b, fl, levs, cnms){
  asgn <- attr(fl, "assign")
  nc <- vapply(cnms, length, 1L)
  nb <- nc * (nl <- vapply(levs, length, 1L)[asgn])
  nbseq <- rep.int(seq_along(nb), nb)
  u <- split(b, nbseq)
  for (i in seq_along(u)){
    u[[i]] <- matrix(u[[i]], ncol = nc[i], byrow = TRUE,
                     dimnames = list(NULL, cnms[[i]]))
  }
  names(u) <- names(cnms)
  u
}
