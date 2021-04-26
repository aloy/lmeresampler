#' @title Plot bootstrap results
#'
#' @description
#' Generate a density plot with a half-eye plot representing the 68% and 95%
#' percentile intervals from an \code{lmeresamp} object.
#'
#' @param x The lmeresamp object to plot.
#' @param var The estimated parameter to plot, as a string or column number.
#' @param ... not used
#'
#' @rdname plot
#' @export 
#' @method plot lmeresamp
#' @importFrom ggplot2 ggplot labs aes
#' @importFrom ggdist stat_halfeye
plot.lmeresamp <- function(x, var, ...){
  
  # set default
  if(missing(var)){
    var <- 1
  }
  
  x$replicates <- as.data.frame(x$replicates)
  
  to_plot <- unlist(x$replicates[var])
  
  ggplot2::ggplot(x$replicates, ggplot2::aes(x = to_plot)) + 
    ggdist::stat_halfeye(fill = "cadetblue", alpha = 0.5) +
    ggplot2::labs(title = paste("density plot of bootstrap estimates for", var), x = var) 
}
