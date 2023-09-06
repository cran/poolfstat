#' S4 class to represent fstats results obtained with computeFstats.
#'
#' @slot f2.values  A data frame with npop(npop-1)/2 rows and 1 (or 3 if blockjackknife is TRUE) columns containing estimates of the f2-statistics over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.)
#' @slot fst.values  A data frame with npop(npop-1)/2 rows and 1 (or 3 if blockjackknife is TRUE) columns containing estimates of the scaled f2.values (same as obtained with compute.pairwiseFST with method="Identity") over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.). The F2 scaling factor is equal to 1-Q2 (where Q2 is the AIS probability between the two populations)
#' @slot f3.values  A data frame with npops(npops-1)(npops-2)/2 rows and 1 (or 4 if blockjackknife is TRUE) columns containing estimates of the f3-statistics over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.) and Z-score measuring the deviation of the f3-statistics from 0 in units of s.e.
#' @slot f3star.values  A data frame with npops(npops-1)(npops-2)/2 rows and 1 (or 4 if blockjackknife is TRUE) columns containing estimates of the scaled f3-statistics over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.) and Z-score measuring the deviation of the f3-statistics from 0 in units of s.e. The F3 scaling factor is equal to 1-Q1 (where Q1 is the AIS probability within the target population, i.e., population C for F3(C;A,B))
#' @slot f4.values  A data frame with npops(npops-1)(npops-2)(npops-3)/8 rows and 1 (or 4 if blockjackknife is TRUE) columns containing estimates of the f4-statistics over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.) and Z-score measuring the deviation of the f4-statistics from 0 in units of s.e.
#' @slot Dstat.values A data frame with npops(npops-1)(npops-2)(npops-3)/8 rows and 1 (or 4 if blockjackknife is TRUE) columns containing estimates of the D-statistics (scaled f4-statistics) over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.) and Z-score measuring the deviation of the f3-statistics from 0 in units of s.e. For a given quadruplet (A,B;C,D), the parameter D corresponds to F4(A,B;C,D) scaled by (1-Q2(A,B))*(1-Q2(C,D)) where Q2(X,Y) is the is the AIS probability between the X and Y populations.
#' @slot F2.bjack.samples If blockjackknife=TRUE and options return.F2.blockjackknife.samples is actived in compute.fstats, an array of dimension (npop x npop x nblocks) in an admixtools2 compatible format
#' @slot comparisons A list containing matrices with population names associated to the different test comparisons (e.g., the "F2" elements of the list is a npop(npop-1)/2 rows x 2 columns with each row containing the name of the two populations compared)
#' @slot Q.matrix The estimated error covariance matrix for all the F2 and F3 estimates (required by graph fitting functions to compute graph scores)
#' @slot heterozygosities  A data frame with npop rows and 1 (or 3 if blockjackknife is TRUE) columns containing estimates of the within population heterozygosities (1-Q1) over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.)
#' @slot divergence  A data frame with npop(npop-1)/2 rows and 1 (or 3 if blockjackknife is TRUE) column(s) containing estimates of each population pairwise (absolute) divergence (1-Q2) over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.). This statistic is related to dXY (a.k.a. PiXY) but it is computed on the ascertained SNPs that were included in the original pooldata or countdata objects.
#' @slot pairwise.fst A npop x npop (symmetric) matrix containing the pairwise-population Fst estimates (same as in the fst.values object) that may directly be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).
#' @slot pairwise.div A npop x npop (symmetric) matrix containing the pairwise-population divergence (1-Q2) estimates (same as in the fst.values object) that may directly be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).
#' @slot blockjacknife A logical indicating whether block-jackknife estimates of standard errors are available (TRUE) or not (FALSE)
#' @seealso To generate pairwise object, see \code{\link{compute.pairwiseFST}}
fstats<-setClass(Class = "fstats",
                 slots=c(f2.values="data.frame",fst.values="data.frame",f3.values = "data.frame",f3star.values = "data.frame",f4.values="data.frame",Dstat.values="data.frame",comparisons="list",Q.matrix="matrix",heterozygosities="data.frame",divergence="data.frame",pairwise.fst="matrix",pairwise.div="matrix",F2.bjack.samples="array",blockjacknife="logical")
)

#' Check fstats objects
#' @param x The name of the object to be tested
#' @export
is.fstats <- function (x) 
{
  res <- (is(x,"fstats") & validObject(x))
  return(res)
}

#' Show fstats object
#' @param object Object of class fstats
setMethod("show","fstats",
          function ( object ){
            cat ( " * * * fstats Object * * * \n" )
            cat("Example of useful visualization functions are plot.fstats")
          }
)

#' plot fstats object
#' @param x Object of class fstats
#' @param y dummy argument
#' @param ... Other arguments to be passed to plot_fstats
#' @seealso see \code{\link{plot_fstats}} for details on plot_fstats arguments
setMethod("plot","fstats",
          function ( x,y, ... ){
            plot_fstats(x,...)
          }
)

#' S4 class to represent a pairwise Fst results obtained with the compute.pairwiseFST
#'
#' @slot values  A data frame with npop*(npop-1)/2 rows and 3 (or 7 if blockjackknife is TRUE) columns containing for both the Fst and Q2, estimates over all the SNPs and if blockjackknife=TRUE, the estimated block-jackknife and standard error (s.e.). The seventh (or third if blockjackknife=FALSE) column gives the number of SNPs.
#' @slot PairwiseFSTmatrix A npxnp matrix containing the pairwise FST estimates
#' @slot PairwiseSnpFST A matrix (nsnp rows and npops columns) with read count data for the reference allele
#' @slot PairwiseSnpQ1 A matrix (nsnp rows and npops columns) with overall read coverage
#' @slot PairwiseSnpQ2 A matrix (nsnp rows and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, allele 1 and allele 2
#' @slot blockjacknife A logical indicating whether block-jackknife estimates of standard errors are available (TRUE) or not (FALSE)
#' @seealso To generate pairwise object, see \code{\link{compute.pairwiseFST}}
pairwisefst<-setClass(Class = "pairwisefst",
                     slots=c(values="data.frame",PairwiseFSTmatrix = "matrix",PairwiseSnpFST="matrix",
                                    PairwiseSnpQ1 = "matrix",PairwiseSnpQ2 = "matrix",blockjacknife="logical")
                   )

#' Check pairwisefst objects
#' @param x The name (or a path) of the pairwisefst object
#' @export
is.pairwisefst <- function (x) 
{
  res <- (is(x,"pairwisefst") & validObject(x))
  return(res)
}

#' Show pairwisefst object
#' @param object Object of class pairwisefst
setMethod("show","pairwisefst",
          function ( object ){
            cat ( " * * * pairwisefst Object * * * \n" )
            cat("Example of useful visualization functions are heatmap, plot.fstats")
          }
)

#' Show pairwisefst object
#' @param x Object of class pairwisefst
#' @param Rowv determines if and how the row dendrogram should be computed and reordered. Either a dendrogram or a vector of values used to reorder the row dendrogram or NA to suppress any row dendrogram (and reordering) or by default, NULL, see ‘Details’ below.
#' @param Colv determines if and how the column dendrogram should be reordered. Has the same options as the Rowv argument above and additionally when x is a square matrix, Colv = "Rowv" means that columns should be treated identically to the rows (and so if there is to be no row dendrogram there will not be a column one either).
#' @param distfun function used to compute the distance (dissimilarity) between both rows and columns. Defaults to as.dist.
#' @param hclustfun function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to hclust. Should take as argument a result of distfun and return an object to which as.dendrogram can be applied.
#' @param reorderfun function(d, w) of dendrogram and weights for reordering the row and column dendrograms. The default uses reorder.dendrogram.
#' @param add.expr expression that will be evaluated after the call to image. Can be used to add components to the plot.
#' @param symm logical indicating if x should be treated symmetrically; can only be true when x is a square matrix.
#' @param revC logical indicating if the column order should be reversed for plotting, such that e.g., for the symmetric case, the symmetry axis is as usual.
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "row" if symm false, and "none" otherwise.
#' @param na.rm logical indicating whether NA's should be removed.
#' @param margins numeric vector of length 2 containing the margins (see par(mar = *)) for column and row names, respectively.
#' @param ColSideColors (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#' @param RowSideColors (optional) character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
#' @param cexRow,cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
#' @param labRow,labCol character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
#' @param main,xlab,ylab main, x- and y-axis titles; defaults to none.
#' @param keep.dendro logical indicating if the dendrogram(s) should be kept as part of the result (when Rowv and/or Colv are not NA).
#' @param verbose logical indicating if information should be printed.
#' @param ... additional arguments passed on to image, e.g., col specifying the colors.
setMethod("heatmap","pairwisefst",
          function ( x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
                     distfun = as.dist, hclustfun = hclust, reorderfun = function(d,w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv,"Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nrow(x@PairwiseFSTmatrix)), cexCol = 0.2 + 1/log10(ncol(x@PairwiseFSTmatrix)), labRow = NULL,labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, verbose = getOption("verbose"), ... ){
             heatmap(x@PairwiseFSTmatrix,Rowv = Rowv, Colv = Colv, 
                     distfun = distfun, hclustfun = hclustfun, reorderfun = reorderfun, add.expr, symm = symm, revC = revC, scale = scale, na.rm = na.rm,margins = margins, ColSideColors, RowSideColors, cexRow = cexRow, cexCol = cexCol, labRow = labRow,labCol = labCol, main = main, xlab = xlab, ylab = ylab, keep.dendro = keep.dendro, verbose =verbose, ...)
          }
)

#' plot pairwisefst object
#' @param x Object of class pairwisefst
#' @param y dummy argument
#' @param ... Some arguments to be passed to plot_fstats
#' @seealso see \code{\link{plot_fstats}} for details on plot_fstats arguments
setMethod("plot","pairwisefst",
          function ( x,y, ... ){
            plot_fstats(x,...)
          }
)
