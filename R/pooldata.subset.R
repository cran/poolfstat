#' Create a subset of the pooldata object that contains Pool-Seq data
#' @param pooldata A pooldata object containing Pool-Seq information
#' @param pool.index Indexes of the pools (at least two), that should be selected to create the new pooldata object
#' @param min.cov.per.pool Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded
#' @param max.cov.per.pool Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overal read counts for the reference allele over the read coverage)
#' @return A pooldata object with 7 elements:
#' \enumerate{
#' \item "refallele.readcount": a matrix with nsnp rows and npools columns containing read counts for the reference allele (chosen arbitrarily) in each pool
#' \item "readcoverage": a matrix with nsnp rows and npools columns containing read coverage in each pool
#' \item "snp.info": a matrix with nsnp rows and four columns containing respectively the contig (or chromosome) name (1st column) and position (2nd column) of the SNP; the allele in the reference assembly (3rd column); the allele taken as reference in the refallele matrix.readcount matrix (4th column); and the alternative allele (5th column)
#' \item "poolsizes": a vector of length npools containing the haploid pool sizes
#' \item "poolnames": a vector of length npools containing the names of the pools
#' \item "nsnp": a scalar corresponding to the number of SNPs
#' \item "npools": a scalar corresponding to the number of pools
#' }
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}}
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  pooldata.subset=pooldata.subset(pooldata,pool.index=c(1,2))
#' @export
pooldata.subset<-function(pooldata,pool.index=c(1,2),min.cov.per.pool=-1,max.cov.per.pool=1e6,min.maf=-1){
  if(!(is.pooldata(pooldata))) {stop("The data are not formatted as a valid pooldata object..(see the readpooldata(), sync2pooldata or vcf2pooldata functions)")}
  npools=length(pool.index)
  if(npools<2){stop("At least two population indexes should be given")}
  if(sum(!(pool.index %in% 1:pooldata@npools))>0){stop("Population indexes should be in the range define by the number of pools of the pooldata object")}

  data.Y=pooldata@refallele.readcount[,pool.index]
  data.N=pooldata@readcoverage[,pool.index]
  ##filtres sur couverture et maf
  tmp.maf=0.5-abs(0.5-rowSums(data.Y)/rowSums(data.N))
  dum.sel=(rowSums(data.N>=min.cov.per.pool)==npools) & (rowSums(data.N<=max.cov.per.pool)==npools) & (tmp.maf>min.maf)
  data.Y=data.Y[dum.sel,] ; data.N=data.N[dum.sel,] 
  
  res<-new("pooldata")
  res@npools=npools
  res@nsnp=nrow(data.Y)
  res@refallele.readcount=data.Y
  rm(data.Y)
  res@readcoverage=data.N
  rm(data.N)
  res@snp.info=pooldata@snp.info[dum.sel,]
  #rm(snpdet)
  res@poolsizes=pooldata@poolsizes[pool.index]
  res@poolnames=pooldata@poolnames[pool.index]
  
  cat("Data consists of",res@nsnp,"SNPs for",res@npools,"Pools\n")
  return(res)
}