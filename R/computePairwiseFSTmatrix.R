#' Compute pairwise population population FST matrix (and possibly all pairwise SNP-specific FST)
#' @param pooldata A pooldata object containing Pool-Seq information
#' @param method Either "Anova" (default method as described in the manuscript) or "Identity" (relies on an alternative modeling consisting in estimating unbiased Probability of Identity within and across pairs of pools)
#' @param min.cov.per.pool Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons.
#' @param max.cov.per.pool Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons.
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overal read counts for the reference allele over the read coverage)  in the pairwise comparisons.
#' @param output.snp.values If TRUE, provide SNP-specific pairwise FST for each comparisons (may lead to a huge result object if the number of pools and/or SNPs is large)
#' @return A list with 2 (or 5 if output.snp.values=TRUE) elements:
#' \enumerate{
#' \item "PairwiseFSTmatrix": a matrix with npools rows and npools columns containing the pairwise pool FST estimates
#' \item "NbOfSNPs": a matrix with npools rows and npools columns containing the number of SNPs satisfying the filtering criteria in pairs of pools (and within each pool in the diagonal)
#' \item "PairwiseSnpFST" (if output.snp.values=TRUE): a matrix with nsnp rows and (npools*(npools-1))/2 columns containing the SNP-specific FST estimates for each pair of pools
#' #' \item "PairwiseSnpQ1" (if output.snp.values=TRUE): a matrix with nsnp rows and (npools*(npools-1))/2 columns containing the SNP-specific Q1 estimates for each pair of pools
#' #' \item "PairwiseSnpQ2" (if output.snp.values=TRUE): a matrix with nsnp rows and (npools*(npools-1))/2 columns containing the SNP-specific Q2 estimates for each pair of pools
#' }
#' @seealso To generate subset of pooldata object, see \code{\link{pooldata.subset}}
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  PairwiseFST=computePairwiseFSTmatrix(pooldata)
#' @export
computePairwiseFSTmatrix<-function(pooldata,method="Anova",min.cov.per.pool=-1,max.cov.per.pool=1e6,min.maf=-1,output.snp.values=FALSE){
  if(!(is.pooldata(pooldata))) {stop("The data are not formatted as a valid pooldata object..(see the readpooldata(), sync2pooldata or vcf2pooldata functions)")}
  if(pooldata@npools<3){stop("At least 3 pools should be given to compute the matrix")}
  
  FSTmatrix=matrix(0,pooldata@npools,pooldata@npools)
  SNPcountmatrix=matrix(0,pooldata@npools,pooldata@npools)
  colnames(FSTmatrix)=rownames(FSTmatrix)=colnames(SNPcountmatrix)=rownames(SNPcountmatrix)=pooldata@poolnames  
  if(output.snp.values){
    snpFST.all=snpQ1.all=snpQ2.all=matrix(NA,pooldata@nsnp,pooldata@npools*(pooldata@npools-1)/2)
    tmp=combn(pooldata@npools,2)
    colnames(snpFST.all)=colnames(snpQ1.all)=colnames(snpQ2.all)=paste0(pooldata@poolnames[tmp[1,]],"_vs_",pooldata@poolnames[tmp[2,]])
  }

  cnt.pair=0  
  for(p in 1:(pooldata@npools-1)){
    for(q in (p+1):(pooldata@npools)){
      cnt.pair=cnt.pair+1
  #on reprend le code de la fonction pooldata.subset car on a besoin d'info sur les SNPs selectionnes (index) au cas ou output.snp.values=T    
      pool.index=c(p,q)
      data.Y=pooldata@refallele.readcount[,pool.index]
      data.N=pooldata@readcoverage[,pool.index]
      ##filtres sur couverture et maf
      OverallN=rowSums(data.N)
      tmp.maf=0.5-abs(0.5-rowSums(data.Y)/OverallN)
      dum.sel=OverallN>0 & (rowSums(data.N>=min.cov.per.pool)==2) & (rowSums(data.N<=max.cov.per.pool)==2) & (tmp.maf>min.maf)
      data.Y=data.Y[dum.sel,] ; data.N=data.N[dum.sel,] 
      
      res<-new("pooldata")
      res@npools=2
      res@nsnp=nrow(data.Y)
      res@refallele.readcount=data.Y
      rm(data.Y)
      res@readcoverage=data.N
      rm(data.N)
      res@snp.info=pooldata@snp.info[dum.sel,]
      #rm(snpdet)
      res@poolsizes=pooldata@poolsizes[pool.index]
      res@poolnames=pooldata@poolnames[pool.index]
      
      fst.res=computeFST(res,method = method)
      
      FSTmatrix[p,q]=FSTmatrix[q,p]=fst.res$FST
      SNPcountmatrix[p,q]=SNPcountmatrix[q,p]=res@nsnp
      if(output.snp.values){
        snpFST.all[dum.sel,cnt.pair]=fst.res$snp.FST
        snpQ1.all[dum.sel,cnt.pair]=fst.res$snp.Q1
        snpQ2.all[dum.sel,cnt.pair]=fst.res$snp.Q2
       }

      cat("Computations for pairs",pooldata@poolnames[p],"vs.",pooldata@poolnames[q],"done\n")                
    }
  }
  
  ##within pool SNP counts
  tmp.maf=0.5-abs(0.5-pooldata@refallele.readcount/pooldata@readcoverage)
  diag(SNPcountmatrix)=colSums(pooldata@readcoverage>=min.cov.per.pool & pooldata@readcoverage<=max.cov.per.pool & tmp.maf>min.maf,na.rm=T)
  
  
  if(output.snp.values){
    res=list(PairwiseFSTmatrix=FSTmatrix,NbOfSNPs=SNPcountmatrix,PairwiseSnpFST=snpFST.all,PairwiseSnpQ1=snpQ1.all,PairwiseSnpQ2=snpQ2.all)  
  }else{
    res=list(PairwiseFSTmatrix=FSTmatrix,NbOfSNPs=SNPcountmatrix)     
  }
  

  return(res)
}