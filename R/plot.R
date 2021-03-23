#' @title plot
#'
#' @description
#' Plots the bootstrap call.
#'
#' @details
#' This function is given \code{x, var} and uses them to plot the density plot
#' of bootstrap estimates for the var of choice.
#'
#' @param x The lmeresamp object to plot.
#' @param var The estimated parameter to plot.
#' @param ... not used
#'
#' @rdname plot
#' @export 
#' @method plot lmeresamp
plot.lmeresamp <- function(x, var, ...){
  
  var_name <- as.character(var)
  
  ggplot(x$replicates, aes(x = x$replicates$var)) + 
    geom_density(fill = "cadetblue", alpha = 0.5) +
    labs(title = "density plot", x = var_name) 
}
