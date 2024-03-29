% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpntPET.R
\name{lpntPET_tol}
\alias{lpntPET_tol}
\title{Changing the tolerance level}
\usage{
lpntPET_tol(obj, tol)
}
\arguments{
\item{obj}{output of the lp_ntPETabc function}

\item{tol}{vector of tolerance levels to check.}
}
\value{
A list of values from the posterior distributions, each element of the list corresponds
to each value of the vector of tolerance levels given in tol.
}
\description{
This function change the automatic tolerance level of the lpntPET and produce
plots comparing posterior distributions with different choices of tolerance level.
Maximum of 6 possible lines to plot.
}
