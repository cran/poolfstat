#' S4 class to represent a Pool-Seq data set.
#'
#' @slot npools The number of pools
#' @slot nsnp The number of SNPs
#' @slot refallele.readcount A matrix (nsnp rows and npools columns) with read count data for the reference allele
#' @slot readcoverage A matrix (nsnp rows and npools columns) with overall read coverage
#' @slot snp.info A data frame (nsnp rows and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, Reference allele name and Alternate allele name (if available)
#' @slot poolsizes A vector of length npools with the corresponding haploid pool sizes
#' @slot poolnames A vector of length npools with the corresponding haploid pool names
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}}, \code{\link{genobaypass2pooldata}} and \code{\link{genoselestim2pooldata}}
pooldata<-setClass(Class = "pooldata",
	slots=c(refallele.readcount = "matrix",readcoverage="matrix",snp.info="data.frame",
	               poolsizes="numeric",poolnames="character",
	               npools = "numeric",nsnp = "numeric")
)

#' Check pooldata objects
#' @param x The name of the object to be tested
#' @export
is.pooldata <- function (x) 
{
    res <- (is(x,"pooldata") & validObject(x))
    return(res)
}

#' Show pooldata object
#' @param object Object of class pooldata
setMethod("show","pooldata",
             function ( object ){
               cat ( " * * * PoolData Object * * * \n" )
               cat ( " * Number of SNPs   = ",object@nsnp,"\n" )
               cat ( " * Number of Pools  = ",object@npools,"\n" )
               cat ( " * Pool Names       :\n",paste(object@poolnames,collapse="; ") )
               cat ( "\n * * * * * * * * * * * * * * \n" )
             }
)