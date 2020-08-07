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
  combo_stats <- sum(map_dbl(boot_list, ~.x$stats)) #change stats
  RES <- boot_list[[1]]
  RES$replicates <- combo_replicates
  RES$stats <- combo_stats
  return(RES)
}