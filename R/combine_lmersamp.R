#' @title Combine bootstrap results
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


#' @title Combine p-values
#' 
#' @description 
#' Combines bootstrap p-values from processes split for parallelization.
#' 
#' @details 
#' This helper function combines the different summary tables produced by
#' \code{bootstrap_pvals()} when run in parallel to yield unified output
#' and a single summary table.
#' 
#' @param ... two or more summary data frames produced by \code{bootstrap_pvals}.
#' 
#' @rdname combine_pvals
#' @export
#' @importFrom tibble as_tibble
combine_pvals <- function(...) {
  n_extreme <- B <- NULL
  pval_list <- list(...)
  pval_tbl <- purrr::map_dfr(
    pval_list, 
    ~dplyr::mutate(.x$coefficients, n_extreme = (p.value * (.x$B + 1)) - 1, B = .x$B)
    )
  
  combo_coef <- pval_tbl %>% 
    dplyr::group_by(dplyr::across(1:4)) %>% 
    # dplyr::group_by(term, Estimate, `Std. Error`, `t value`) %>% 
    dplyr::summarize(p.value = (sum(n_extreme) + 1) / (sum(B) + 1), .groups = "drop_last") %>%
    dplyr::ungroup()
  
  combo_seeds <- purrr::map(pval_list, ~.x$seed)
  
  combo_r <- sum(map_dbl(pval_list, ~.x$B))
  type <- unique(sapply(pval_list, function(x) x$type))
    
  structure(
    list(coefficients = combo_coef, B = combo_r, seed = combo_seeds, type = type),
    class = "coef_tbl"
  )
}
