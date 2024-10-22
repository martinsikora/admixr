% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/adm.R
\name{llAdmixture}
\alias{llAdmixture}
\title{likelihood function for ADMIXTURE projections}
\usage{
llAdmixture(qM, pM, gt)
}
\arguments{
\item{qM}{The Q matrix of cluster proportions from ADMIXTURE,
dimensions individuals x components}

\item{pM}{The P matrix of ancestral cluster frequencies from ADMIXTURE,
dimensions n SNPs x components}

\item{gt}{genotype vector}
}
\value{
log-likelihood
}
\description{
\code{llAdmixture} calculates the likelihood of the genotypes
given the Q matrix and P matrix.
genotype vector is coded as 0,1,2 and has to match coding from P matrix
Equation (2) from Alexander et al. Genome Res. 2009 19: 1655-1664
}

