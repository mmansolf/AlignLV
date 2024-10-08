% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{align.optim}
\alias{align.optim}
\title{Runs alignment optimizer}
\usage{
align.optim(
  stacked,
  n,
  estimator,
  nstarts = 50,
  ncores = 3,
  hyper.first,
  center.means,
  eps.alignment,
  clf.ignore.quantile,
  verbose
)
}
\arguments{
\item{stacked}{Stacked parameter estimates from \code{\link{stackEstimates}}}

\item{n}{Sample size in each group}

\item{estimator}{See \code{\link{Alignment}} documentation.}

\item{nstarts}{Number of starting values for alignment; default is 10}

\item{ncores}{See \code{\link{Alignment}} documentation.}

\item{hyper.first}{See \code{\link{Alignment}} documentation.}

\item{center.means}{See \code{\link{Alignment}} documentation.}

\item{eps.alignment}{See \code{\link{Alignment}} documentation.}

\item{clf.ignore.quantile}{See \code{\link{Alignment}} documentation.}

\item{verbose}{See \code{\link{Alignment}} documentation.}
}
\value{
A \code{list} of results from multiple runs of the alignment optimizer:
\itemize{
\item \code{mv} Means and variances from each alignment run.
\item \code{parout} A table of outputs from \code{link[stats]{optim}} containing the
function values, convergence information, and resulting estimates of means
and variances from each run.
\item \code{nFailedRuns} The number of runs that failed to complete. An error is
returned if no runs fail.
}
}
\description{
Not generally intended to be used on its own, but exported anyway
for didactic purposes.
}
\details{
See example for \code{\link{Alignment}} for examples
}
