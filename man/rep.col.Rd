% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpntPET.R
\name{rep.col}
\alias{rep.col}
\title{Repeating columns}
\usage{
\method{rep}{col}(x, n)
}
\arguments{
\item{x}{vector to be repeated. The length of this vector will be the number of rows
of the matrix.}

\item{n}{scalar. It indicates the number of times the vector is repeated. It will be
the number of columns of the matrix.}
}
\value{
A matrix with dimension length(x) x n
}
\description{
This function repeats the same vector to create a matrix with n repeated columns.
}
