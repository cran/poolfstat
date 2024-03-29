% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_blockDdenom}
\alias{compute_blockDdenom}
\alias{.compute_blockDdenom}
\title{compute_blockDdenom}
\usage{
.compute_blockDdenom(refcount, totcount, nblocks, block_id, verbose)
}
\arguments{
\item{refcount}{Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele}

\item{totcount}{Matrix of nsnpxnpop with total counts or read coverages}

\item{nblocks}{Integer giving the number of block-jackknife blocs (may be 0 if no block-jackknife)}

\item{block_id}{Integer vector of length nsnps with the (0-indexed) id of the block to which each SNP belongs (-1 for SNPs outside blocks)}

\item{verbose}{Logical (if TRUE progression bar is printed on the terminal)}
}
\value{
Return a matrix with nf4=(npops*(npops-1)/2)*((npops-2)*(npops-3)/2)/2 rows and nblocks+1 columns giving the mean Dstat-denominator (1-Q2ab)(1-Q2cd)
 for all quadruplet configuration and within each block-jackknife sample and over all SNPs (last column)
}
\description{
Compute the denominator of the Dstat for all quadruplet configuration and each block-jackknife block (if any) and overall SNPs (within or outside blocks)
}
\details{
Compute the denominator of the Dstat for all quadruplet configuration and each block-jackknife block (if any) and overall SNPs (within or outside blocks)
}
\examples{
#
}
