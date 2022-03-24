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
#' @importFrom tidyr pivot_longer
plot.lmeresamp <- function(x, var, ...){
  value <- term <- NULL
  if(is.numeric(x$replicates)) {
    ggplot2::ggplot(data = NULL, ggplot2::aes(x = x$replicates)) + 
      ggdist::stat_halfeye(fill = "cadetblue", alpha = 0.5)
  } else{
  
    # set default
    if(missing(var)){
      tidy_reps <- tidyr::pivot_longer(
        x$replicates, 
        cols = dplyr::everything(), 
        names_to = "term", 
        values_to = "value"
      )
      
      ggplot2::ggplot(tidy_reps, ggplot2::aes(x = value, y = term)) + 
        ggdist::stat_halfeye(fill = "cadetblue", alpha = 0.5)
      
    } else{
      x$replicates <- as.data.frame(x$replicates)
      if(is.numeric(var)) var <- colnames(x$replicates)[var]
      if(grepl("[()]", var)) var <- paste0("`", var, "`")
      
      # to_plot <- unlist(x$replicates[var])
      
      ggplot2::ggplot(x$replicates, ggplot2::aes_string(x = var)) + 
        ggdist::stat_halfeye(fill = "cadetblue", alpha = 0.5) +
        ggplot2::labs(
          title = paste("Distribution of", var), 
          x = var,
          y = "density"
        ) 
    }
    
  }
}
