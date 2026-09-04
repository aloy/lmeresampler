.onLoad <- function(libname, pkgname) {
  register_s3_method("lme4", "isGLMM", "glmmTMB", isGLMM.glmmTMB)
}

#' Register an S3 method for a generic owned by an optional dependency
#'
#' @description
#' Registers \code{method} against \code{pkg::generic} for \code{class}
#' once \code{pkg} is loaded, without forcing \code{pkg} (an optional,
#' \code{Suggests}-only dependency) to load just because this package is
#' loaded. See
#' \url{https://r-pkgs.org/dependencies-in-practice.html#sec-dependencies-suggests-conditional}.
#' @keywords internal
#' @noRd
register_s3_method <- function(pkg, generic, class, method) {
  register <- function(...) {
    envir <- asNamespace(pkg)
    if (exists(generic, envir = envir, inherits = FALSE)) {
      registerS3method(generic, class, method, envir = envir)
    }
  }

  if (pkg %in% loadedNamespaces()) {
    register()
  }

  setHook(packageEvent(pkg, "onLoad"), register)
}
