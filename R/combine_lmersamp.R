#' @title combine
#'
#' @description
#' Combines processes split for parallelization.
#'
#' @details
#' This helper function combines the different processes split for 
#' parallelization to yield unified output and bootstrap statistics.
#'
#' @rdname combine
#' @export 
# bootstrap CI method for object of class lmeresamp
combine_lmeresamp <- function(...) {
  boot_list <- list(...)
  combo_replicates <- purrr::map_dfr(boot_list, ~as.data.frame(.x$replicates))
  combo_r <- sum(map_dbl(boot_list, ~.x$R))
  RES <- boot_list[[1]]
  RES$R <- combo_r
  RES$replicates <- combo_replicates
  RES$stats$rep.mean <- colMeans(RES$replicates) # recalculated mean
  RES$stats$se <- unlist(purrr::map(RES$replicates, sd)) # recalculated se
  RES$stats$bias <- RES$stats$rep.mean - RES$stats$observed # recalculated bias
  return(RES)
}

