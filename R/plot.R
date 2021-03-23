#' @title plot
#'
#' @description
#' Plots the bootstrap call.
#'
#' @details
#' This function is given \code{x, var} and uses them to plot the density plot
#' of bootstrap estimates for the var of choice.
#'
#' @param obj The lmeresamp object to plot.
#' @param var The estimated parameter to plot, as a string.
#' @param ... not used
#'
#' @rdname plot
#' @export 
#' @method plot lmeresamp
plot.lmeresamp <- function(obj, var, ...){
  
  obj$replicates <- as.data.frame(obj$replicates)
  
  to_plot <- unlist(obj$replicates[var])
  
  ggplot2::ggplot(obj$replicates, ggplot2::aes(x = to_plot)) + 
    ggplot2::geom_density(fill = "cadetblue", alpha = 0.5) +
    ggplot2::labs(title = "density plot", x = var) 
}
