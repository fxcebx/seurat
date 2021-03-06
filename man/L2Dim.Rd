% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensional_reduction.R
\name{L2Dim}
\alias{L2Dim}
\title{L2-normalization}
\usage{
L2Dim(object, reduction, new.dr = NULL, new.key = NULL)
}
\arguments{
\item{object}{Seurat object}

\item{reduction}{Dimensional reduction to normalize}

\item{new.dr}{name of new dimensional reduction to store
(default is olddr.l2)L2Norm}

\item{new.key}{name of key for new dimensional reduction}
}
\value{
Returns a \code{\link{Seurat}} object
}
\description{
Perform l2 normalization on given dimensional reduction
}
