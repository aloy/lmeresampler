% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_lmersamp.R
\name{combine_lmeresamp}
\alias{combine_lmeresamp}
\title{Combine bootstrap results}
\usage{
combine_lmeresamp(...)
}
\arguments{
\item{...}{two or more objects of class \code{lmeresamp}, to be combined into one.}
}
\description{
Combines bootstrap results from processes split for parallelization.
}
\details{
This helper function combines the different processes split for
parallelization to yield unified output and bootstrap statistics.
}
