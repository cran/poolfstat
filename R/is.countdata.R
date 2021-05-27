#' S4 class to represent a Count data set.
#'
#' @slot npops The number of populations
#' @slot nsnp The number of SNPs
#' @slot refallele.count A matrix (nsnp rows and npops columns) with the allele counts for the reference allele
#' @slot total.count A matrix (nsnp rows and npops columns) with the total number of counts (i.e., twice the number of genotyped individual for diploid species and autosomal markers)
#' @slot snp.info A data frame (nsnp rows and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, Reference allele name and Alternate allele name (if available)
#' @slot popnames A vector of length npops with the corresponding population names
#' @seealso To generate countdata object, see \code{\link{genobaypass2countdata}} and \code{\link{genotreemix2countdata}}
countdata<-setClass(Class = "countdata",
                   slots = c(refallele.count = "matrix",total.count="matrix",snp.info="data.frame",
                                  popnames="character",npops = "numeric",nsnp = "numeric")
)

#' Check countdata objects
#' @param x The name of the object to be tested
#' @export
is.countdata <- function (x) 
{
  res <- (is(x,"countdata") & validObject(x))
  return(res)
}

#' Show countdata object
#' @param object Object of class countdata
setMethod("show","countdata",
          function ( object ){
            cat ( " * * * Countdata Object * * * \n" )
            cat ( " * Number of SNPs  = ",object@nsnp,"\n" )
            cat ( " * Number of Pops  = ",object@npops,"\n" )
            cat ( " * Pop Names       :\n",paste(object@popnames,collapse="; ") )
            cat ( "\n * * * * * * * * * * * * * * \n" )
          }
)