% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/fstats.R
\name{doF3Test}
\alias{doF3Test}
\title{F3 test}
\usage{
doF3Test(p, idxM, blockIdx)
}
\arguments{
\item{p}{allele counts array, dimensions SNPs x individuals x 2}

\item{idxM}{matrix of dimensions 1 x 3 with population labels for test
configuration f3(test;ref1,ref2)}

\item{blockIdx}{index vector of chromosome blocks for jackknife}
}
\value{
a data frame with f3 results
\itemize{
  \item p1, p2, p3  .. test configuration
  \item f3 .. estimate
  \item se, Z .. standard error and Z score from block jackknife
  \item nSites .. number of informative SNPs used in the test
}
}
\description{
Takes an array of allele counts input and does f3 test as in
Patterson et al. 2012 on specified input populations; corresponds
to "f3 outgroup" statistic if 'test' is an outgroup population
}

