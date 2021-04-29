#' @title combine
#'
#' @description
#' Combines bootstrap results from processes split for parallelization.
#'
#' @details
#' This helper function combines the different processes split for
#' parallelization to yield unified output and bootstrap statistics.
#'
#' @param ... two or more objects of class \code{lmeresamp}, to be combined into one.
#'
#' @rdname combine
#' @export 
#' @importFrom purrr map_dfr map map_dbl
combine_lmeresamp <- function(...) {
  boot_list <- list(...)
  combo_replicates <- purrr::map_dfr(boot_list, ~as.data.frame(.x$replicates))
  combo_r <- sum(map_dbl(boot_list, ~.x$B))
  RES <- boot_list[[1]]
  RES$B <- combo_r
  RES$replicates <- combo_replicates
  RES$stats$rep.mean <- colMeans(RES$replicates) # recalculated mean
  RES$stats$se <- unlist(purrr::map(RES$replicates, sd)) # recalculated se
  RES$stats$bias <- RES$stats$rep.mean - RES$stats$observed # recalculated bias
  return(RES)
}


