% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{gen_true_RHO}
\alias{gen_true_RHO}
\title{gen_true_RHO}
\usage{
gen_true_RHO(wRHO = 1, NN, QQ, RHO = NULL)
}
\arguments{
\item{wRHO}{An integer taking values 1, 2, or 3 for one of the three
scenarios. For \code{wRHO} equal to 2 or 3, \code{QQ} needs to be 
set to 3 cell types.}

\item{NN}{Positive integer for sample size.}

\item{QQ}{Positive integer for number of cell types.}

\item{RHO}{Default to \code{NULL} leads to simulating cell 
type proportions. If a matrix of proportions is supplied, 
this function will append row/column names.}
}
\description{
Simulates cell type proportions under three scenarios.
}
