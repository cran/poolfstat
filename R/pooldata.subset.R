#' Create a subset of the pooldata object that contains Pool-Seq data as a function of pool and/or SNP indexes
#' @param pooldata A pooldata object containing Pool-Seq information
#' @param pool.index Indexes of the pools (at least two), that should be selected to create the new pooldata object (default=all the pools)
#' @param snp.index Indexes of the SNPs (at least two), that should be selected to create the new pooldata object (default=all the SNPs)
#' @param min.cov.per.pool Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded
#' @param max.cov.per.pool Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded
#' @param cov.qthres.per.pool A two-elements vector containing the minimal (qmin) and maximal (qmax) quantile coverage thresholds applied to each pools (0<=qmin<qmax<=1). See details below
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio over all read counts for the reference allele over the read coverage)
#' @param return.snp.idx If TRUE, the row.names of the snp.info slot of the returned pooldata object are named as "rsx" where x is the index of SNP in the initial pooldata object (default=FALSE)
#' @param verbose If TRUE return some information
#' @details This function allows subsetting a pooldata object by selecting only some pools and/or some SNPs (e.g., based on their position on the genome). Additional filtering steps on SNPs can be carried out on the resulting subset to discard SNP with low polymorphism or poorly or too highly covered. In addition, coverage criteria can be applied on a per-pool basis with the cov.qthres.per.pool argument. 'more specific SNP selection based on their positions on the genome or their characteristics. For instance if qmax=0.95, a position is discarded if in a given pool it has a number of reads higher than the 95-th percentile of the empirical coverage distribution in this same pool (defined over the SNPs selected by snp.index). Similarly, if qmax=0.05, a position is discarded if in a given pool it has a number of reads lower than the 5-th percentile of the empirical coverage distribution in this same pool. This mode of selection may be more relevant when considering pools with heterogeneous read coverages. 
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
#'  subset.by.pools=pooldata.subset(pooldata,pool.index=c(1,2))
#'  subset.by.snps=pooldata.subset(pooldata,snp.index=10:100)
#'  subset.by.pools.and.snps=pooldata.subset(pooldata,pool.index=c(1,2),snp.index=10:100)
#'  subset.by.pools.qcov.thr=pooldata.subset(pooldata,pool.index=1:8,cov.qthres.per.pool=c(0.05,0.95)) 
#' @export
pooldata.subset<-function(pooldata,pool.index=1:pooldata@npools,snp.index=1:pooldata@nsnp,min.cov.per.pool=-1,max.cov.per.pool=1e6,min.maf=-1,cov.qthres.per.pool=c(0,1),return.snp.idx=FALSE,verbose=TRUE){
  if(!(is.pooldata(pooldata))) {stop("The data are not formatted as a valid pooldata object (see the vcf2pooldata() or similar functions\n)")}
  if(cov.qthres.per.pool[1]<0 | cov.qthres.per.pool[2]>1){stop("cov.qthres.per.pool must range between 0 and 1\n")}
  if(cov.qthres.per.pool[1]>cov.qthres.per.pool[2]){stop("The two elements in cov.qthres.per.pool must be increasing order\n")}  
  npools=length(pool.index)
  if(npools<2){stop("At least two population indexes should be given\n")}
  if(max(pool.index)>pooldata@npools){stop("Pool indexes should be in the range defined by the number of pools in the pooldata object")}
  if(max(snp.index)>pooldata@nsnp){stop("SNP indexes should be in the range defined by the number of SNPs in the pooldata object")}
  data.Y=pooldata@refallele.readcount[snp.index,pool.index]
  data.N=pooldata@readcoverage[snp.index,pool.index]
  data.D=pooldata@snp.info[snp.index,]
  if(return.snp.idx){rownames(data.D)=paste0("rs",snp.index)} #useful for computePairwiseFSTmatrix
  ##filtres sur couverture et maf
  OverallN=rowSums(data.N)
  tmp.maf=0.5-abs(0.5-rowSums(data.Y)/OverallN)
  dum.sel=OverallN>0 & (rowSums(data.N>=min.cov.per.pool)==npools) & (rowSums(data.N<=max.cov.per.pool)==npools) & (tmp.maf>min.maf) 
  data.Y=data.Y[dum.sel,] ;  data.N=data.N[dum.sel,] ;  data.D=data.D[dum.sel,]
  ##coverage per pools if needed
  if(cov.qthres.per.pool[1]>0 | cov.qthres.per.pool[2]<1){
    min.depth=apply(data.N,2,quantile,probs=cov.qthres.per.pool[1])    
    max.depth=apply(data.N,2,quantile,probs=cov.qthres.per.pool[2])
    index.pos.sel=rowSums(data.N>(rep(1,nrow(data.N))%*%t(max.depth)))==0 & rowSums(data.N<(rep(1,nrow(data.N))%*%t(min.depth)))==0
    if(sum(index.pos.sel)==1){
      data.Y=matrix(data.Y[index.pos.sel,],nrow=1)
      data.N=matrix(data.N[index.pos.sel,],nrow=1)
    }else{
      data.Y=data.Y[index.pos.sel,]
      data.N=data.N[index.pos.sel,]  
    }
    data.D=data.D[index.pos.sel,] 
  }  
  res<-new("pooldata")
  res@npools=npools
  res@nsnp=nrow(data.Y)
  res@refallele.readcount=data.Y
  rm(data.Y)
  res@readcoverage=data.N
  rm(data.N)
  res@snp.info=data.D
  res@poolsizes=pooldata@poolsizes[pool.index]
  res@poolnames=pooldata@poolnames[pool.index]
  
  if(verbose){cat("Data consists of",res@nsnp,"SNPs for",res@npools,"Pools\n")}
  return(res)
}