#' Create a subset of a countdata object that contains count data as a function of pop or SNP indexes
#' @param countdata A countdata object containing Allele count information
#' @param pop.index Indexes of the pops (at least two), that should be selected to create the new pooldata object (default=all the pops)
#' @param snp.index Indexes of the SNPs (at least two), that should be selected to create the new pooldata object (default=all the SNPs)
#' @param min.indgeno.per.pop  Minimal number of overall counts required in each population. If at least one pop is not genotyped for at least min.indgeno.per.pop (haploid) individual, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overall counts for the reference allele over the overall number of (haploid) individual genotyped)
#' @param return.snp.idx If TRUE, the row.names of the snp.info slot of the returned pooldata object are named as "rsx" where x is the index of SNP in the initial pooldata object (default=FALSE)
#' @param verbose If TRUE return some information
#' @details This function allows subsetting a pooldata object by selecting only some pools and/or some SNPs (e.g., based on their position on the genome). Additional filtering steps on SNPs can be carried out on the resulting subset to discard SNP with low polymorphism or poorly or too highly covered. In addition, coverage criteria can be applied on a per-pool basis with the cov.qthres.per.pool argument. 'more specific SNP selection based on their positions on the genome or their characteristics. For instance if qmax=0.95, a position is discarded if in a given pool it has a number of reads higher than the 95-th percentile of the empirical coverage distribution in this same pool (defined over the SNPs selected by snp.index). Similarly, if qmax=0.05, a position is discarded if in a given pool it has a number of reads lower than the 5-th percentile of the empirical coverage distribution in this same pool. This mode of selection may be more relevant when considering pools with heterogeneous read coverages. 
#' @return A countdata object with 6 elements:
#' \enumerate{
#' \item "refallele.count": a matrix (nsnp rows and npops columns) with the allele counts for the reference allele
#' \item "total.count": a matrix (nsnp rows and npops columns) with the total number of counts (i.e., twice the number of genotyped individual for diploid species and autosomal markers)
#' \item "snp.info": a matrix with nsnp rows and four columns containing respectively the contig (or chromosome) name (1st column) and position (2nd column) of the SNP; the allele taken as reference in the refallele.count matrix (3rd column); and the alternative allele (4th column)
#' \item "popnames": a vector of length npops containing the names of the pops
#' \item "nsnp": a scalar corresponding to the number of SNPs
#' \item "npops": a scalar corresponding to the number of populations
#' }
#' @seealso To generate countdata object, see \code{\link{genobaypass2countdata}}, \code{\link{genotreemix2countdata}}
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  pooldata2genobaypass(pooldata=pooldata,writing.dir=tempdir())
#'  ##NOTE: This example is just for the sake of illustration as it amounts to
#'  ##interpret read count as allele count which must not be done in practice!
#'  countdata=genobaypass2countdata(genobaypass.file=paste0(tempdir(),"/genobaypass")) 
#'  subset.by.snps=countdata.subset(countdata,snp.index=10:100)
#'  subset.by.pops.and.snps=countdata.subset(countdata,pop.index=c(1,2),snp.index=10:100)
#' @export
countdata.subset<-function(countdata,pop.index=1:countdata@npops,snp.index=1:countdata@nsnp,min.indgeno.per.pop=-1,min.maf=-1,return.snp.idx=FALSE,verbose=TRUE){
  if(!(is.countdata(countdata))) {stop("The data are not formatted as a valid countdata object (see the genobaypass2countdata() or genotreemix2countdata functions\n)")}
  npops=length(pop.index)
  if(npops<2){stop("At least two population indexes should be given\n")}
  if(max(pop.index)>countdata@npops){stop("Pop indexes should be in the range defined by the number of pops in the countdata object")}
  if(max(snp.index)>countdata@nsnp){stop("SNP indexes should be in the range defined by the number of SNPs in the countdata object")}
  data.Y=countdata@refallele.count[snp.index,pop.index]
  data.N=countdata@total.count[snp.index,pop.index]
  data.D=countdata@snp.info[snp.index,]
  if(return.snp.idx){rownames(data.D)=paste0("rs",snp.index)} #useful for computePairwiseFSTmatrix
  ##filtres sur couverture et maf
  OverallN=rowSums(data.N)
  tmp.maf=0.5-abs(0.5-rowSums(data.Y)/OverallN)
  dum.sel=OverallN>0 & (rowSums(data.N>=min.indgeno.per.pop)==npops) & (tmp.maf>min.maf) 
  res<-new("countdata")
  res@npops=npops
  res@nsnp=sum(dum.sel)
  res@refallele.count=data.Y[dum.sel,]
  rm(data.Y) ; gc()
  res@total.count=data.N[dum.sel,]
  rm(data.N) ; gc()
  res@snp.info=data.D[dum.sel,]
  res@popnames=countdata@popnames[pop.index]
  
  if(verbose){cat("Data consists of",res@nsnp,"SNPs for",res@npops,"Pops\n")}
  return(res)
}