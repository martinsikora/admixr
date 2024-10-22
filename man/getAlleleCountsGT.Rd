% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/utils.R
\name{getAlleleCountsGT}
\alias{getAlleleCountsGT}
\title{Get allele counts from a genotype matrix}
\usage{
getAlleleCountsGT(gt, popIdx, isHaploid = FALSE)
}
\arguments{
\item{gt}{A matrix of genotypes, dimensions SNPs x individuals}

\item{popIdx}{A vector specifying the population label for each individual (column) of the genotype matrix}

\item{isHaploid}{Optional flag if genotype matrix represents haploid data coded as homozygous diploid, i.e. only 0 or 2 entries}
}
\value{
An array with dimensions (SNPs x Populations x 2)
}
\description{
Takes a  genotype matrix coded as 0,1,2 as input and returns an
array with the counts of the two alleles
for each SNP in each population
}

