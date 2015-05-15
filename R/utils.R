#' Tabulate occurences in a vector
#'
#' @param x A vector of values
#' @return A named vector with occurences of the elements of x, sorted along x
#' @examples
#' fastTable(c(1:5, 2))
#' fastTable(c("D", "A", "C", "C"))
fastTable <- function(x){
    uniqueX <- sort(unique(x))
    idx <- match(x, uniqueX)
    res <- tabulate(idx)
    names(res) <- uniqueX
    return(res)
}


#' Generate index vector for chromosome blocks
#'
#' Takes vectors of chromosome IDs and positions as inputs, together with block size
#' returns an index vector corresponding to segments of desired size
#' @param chromosomes A vector of chromosome IDs
#' @param positions A vector of chromosomal positions
#' @param blockSize size of the blocks
#' @return A vector with the block index for each input position
getChromosomeBlocks <- function(chromosomes, positions, blockSize){
    chromosomes <- factor(chromosomes, levels = unique(chromosomes)) ## to keep chromosomes in input order for following 'split' command
    d <- split(positions, chromosomes)

    splitIdx <- lapply(d, function(x) {
        x %/% blockSize
    }) ## integer index of blocks per chromosome
    d1 <- cumsum(sapply(splitIdx, max) + 1) ## max index value per chromosome
    d1 <- d1 - min(d1) ## offset to add for each chromosome (i.e. )
    res <- unlist(mapply("+", splitIdx, d1))
    return(res)
}


#' Get allele counts from a genotype matrix
#'
#' Takes a  genotype matrix coded as 0,1,2 as input and returns an
#' array with the counts of the two alleles
#' for each SNP in each population
#'
#' @param gt A matrix of genotypes, dimensions SNPs x individuals
#' @param popIdx A vector specifying the population label for each individual (column) of the genotype matrix
#' @param isHaploid Optional flag if genotype matrix represents haploid data coded as homozygous diploid, i.e. only 0 or 2 entries
#' @return An array with dimensions (SNPs x Populations x 2)
getAlleleCountsGT <- function(gt, popIdx, isHaploid = FALSE){

    ## helpers
    popIds <- unique(popIdx)
    nPops <- length(popIds)
    nMarkers <- nrow(gt)

    ## number of nonmissing genotypes
    gtM <- matrix(as.integer(!is.na(gt)), ncol = ncol(gt))
    nObs <- t(rowsum(t(gtM), popIdx, reorder = FALSE)) * as.integer(2)

    ## number of alternative alleles per pop & snp
    nAlt <- t(rowsum(t(gt), popIdx, reorder = FALSE, na.rm = TRUE))
    nAlt[nObs == 0] <- NA ## set alt allele count for SNPs without data to NA

    ## build array of counts for marker * population * allele
    counts <- array(NA, dim = c(nMarkers, nPops, 2), dimnames = list(marker = rownames(gt), population = popIds, allele = c("allele1", "allele2")))
    counts[,,2] <- nAlt
    counts[,,1] <- nObs - counts[,,2] ## nObs == 0 SNPs will be set to NA here, since nAlt for those is NA
    if(isHaploid){
        counts <- counts %/% as.integer(2) ## integer division to avoid conversion to numeric
    }
    return(counts)
}

