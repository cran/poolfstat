#' An S4 class to represent a Pool-Seq data set.
#'
#' @slot npools The number of pools
#' @slot nsnp The number of SNPs
#' @slot refallele.readcount A matrix (nsnp rows and npools columns) with read count data for the reference allele
#' @slot readcoverage A matrix (nsnp rows and Npools columns) with overall read coverage
#' @slot snp.info A matrix (nsnp rows and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, the reference allele and the alternative allele
#' @slot poolsizes A vector of length npools with the corresponding haploid pool sizes
#' @slot poolnames A vector of length npools with the corresponding haploid pool names
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}}, \code{\link{genobaypass2pooldata}} and \code{\link{genoselestim2pooldata}}
pooldata<-setClass(Class = "pooldata",
	representation(refallele.readcount = "matrix",readcoverage="matrix",snp.info="matrix",
	               poolsizes="numeric",poolnames="character",
	               npools = "numeric",nsnp = "numeric")
)

#' Check pooldata objects
#' @param x The name (or a path) of the Popoolation sync file (might be in compressed format)
is.pooldata <- function (x) 
{
    res <- (is(x,"pooldata") & validObject(x))
    return(res)
}