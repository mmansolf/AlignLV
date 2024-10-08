% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{getEstimates.lavaan}
\alias{getEstimates.lavaan}
\title{Prepare \code{lavaan} estimates for alignment}
\usage{
getEstimates.lavaan(fit, SE = TRUE)
}
\arguments{
\item{fit}{A \code{lavaan} object compatible with \code{\link{Alignment}}}

\item{SE}{logical; whether to also obtain standard errors.}
}
\value{
A \code{list} of estimates in a format amenable to subsequent alignment
}
\description{
Not generally intended to be used on its own, but exported anyway for
didactic purposes.
}
\details{
See example for \code{\link{Alignment}} for examples

This program was designed based on the published work of Asparouhov &
Muthen, and was not intended to match Mplus exactly, and may not.
}
