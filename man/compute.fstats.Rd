% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.fstats.R
\name{compute.fstats}
\alias{compute.fstats}
\title{Estimate the F-statistics (F2, F3, F3star, F4, Dstat) and within and across population diversity}
\usage{
compute.fstats(
  x,
  nsnp.per.bjack.block = 0,
  computeDstat = FALSE,
  computeF3 = TRUE,
  computeF4 = TRUE,
  output.pairwise.fst = TRUE,
  output.pairwise.div = TRUE,
  computeQmat = TRUE,
  return.F2.blockjackknife.samples = FALSE,
  return.F4.blockjackknife.samples = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{x}{A pooldata object containing Pool-Seq information or a countdata object containing allele count information}

\item{nsnp.per.bjack.block}{Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling)}

\item{computeDstat}{If TRUE compute Dstatistics (i.e. scaled F4). This may add some non negligible computation time if the number of population is large (n>15)}

\item{computeF3}{If TRUE (default) compute all F3 and all F3star (i.e. scaled F3).}

\item{computeF4}{If TRUE (default) compute all F4.}

\item{output.pairwise.fst}{If TRUE (default), output the npopxnpop matrix of pairwise-population Fst estimates (corresponding to the "Identity" method implemented in \code{\link{compute.pairwiseFST}}) in the pairwise.fst slot of the fstats output object (see help(fstats) for details) that may be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).}

\item{output.pairwise.div}{If TRUE (default), output the npopxnpop matrix of pairwise-population divergence (1-Q2) estimates  in the pairwise.div slot of the fstats output object (see help(fstats) for details) that may be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).}

\item{computeQmat}{If TRUE, compute the error covariance matrix between all F3 and F2 statistics (needed for admixture graph construction). This matrix may be very large if the number of pops is large. It is recommended to estimate it on a reduced sample of pops.}

\item{return.F2.blockjackknife.samples}{If TRUE (and nsnp.per.bjack.block>0) return an array of dimension (npopxnpopxnblocks) in an admixtools2 compatible format}

\item{return.F4.blockjackknife.samples}{Deprecated options (since v. 2.2.0)}

\item{verbose}{If TRUE extra information is printed on the terminal}
}
\value{
An object of class fstats (see help(fstats) for details)
}
\description{
Estimate the F-statistics (F2, F3, F3star, F4, Dstat) and within and across population diversity
}
\details{
The function estimates for the n populations (or pools) represented in the input object x:
\enumerate{
\item The F2 statistics for all the \eqn{n(n-1)/2} pairs of populations (or pools) and their scaled version (equivalent, but faster, than Fst estimated with \code{\link{compute.pairwiseFST}} when method="Identity")
\item If n>2, The F3 statistics for all the \eqn{npools(npools-1)(npools-2)/2} possible triplets of populations (or pools) and their scaled version (named F3star after Patterson et al., 2012)
\item If n>3, The F4 statistics and the D-statistics (a scaled version of the F4) for all the \eqn{npools(npools-1)(npools-2)*(npools-3)/8} possible quadruplets of populations
\item The estimated within population heterozygosities (=1-Q1)
\item The estimated divergence for each pair of populations (=1-Q2) 
}
}
\examples{
 make.example.files(writing.dir=tempdir())
 pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
 res.fstats=compute.fstats(pooldata)
}
\seealso{
To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
}
