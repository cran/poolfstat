#' PCA of a pooldata or countdata object using a random allele approach
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param scale If FALSE the random allele data matrix is not scaled (default=TRUE)
#' @param return.snploadings If TRUE return the SNP loadings (may be large)
#' @param plot.pcs A vector with two-elements giving the two PCs to plot. If NULL, no plotting is done.
#' @param ... graphical parameters (see \code{\link{plot}} function)
#' @details PCA is performed by singular-value decomposition (SVD) of a npop (or npools) x nsnp matrix of a single randomly sampled allele (i.e. or read for pooldata object) for each SNP and for each population (inspired by Skoglund and Jakobsson, 2011, https://doi.org/10.1073/pnas.1108181108). Although this approach leads to information loss, it allows to efficiently account for unequal sample size (and read coverages for pool-seq data) and have little impact on the resulting representation when the number of SNPs is large. Note also that the implemented approach is similar to that implemented in the PCA_MDS module of the software ANGSD by Korneliussen et al. (2014) (see http://www.popgen.dk/angsd/index.php/PCA_MDS).
#' @return An object of class fstats (see help(fstats) for details)
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata<-popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.pca<-randomallele.pca(pooldata)
#' @export
randomallele.pca<-function(x,scale=TRUE,return.snploadings=FALSE,plot.pcs=c(1,2),...){
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")} 
  if(x@nsnp<1){stop("No SNP available")}
  if(is.pooldata(x)){popnames=x@poolnames}else{popnames=x@popnames}
  npop=length(popnames)
  if(!is.null(plot.pcs)){
  if(length(plot.pcs)!=2 | min(plot.pcs)<1 | max(plot.pcs)>min(npop-1,x@nsnp-1)){
    stop("Invalid plot.pcs argument (should either be NULL or a two-elements vector with the two PC indexes to be plotted (which must be integers >0 and <npops)")}
  }
  #allele frequencies
  if(is.pooldata(x)){
    x.svd=x@refallele.readcount/x@readcoverage
  }else{
    x.svd=x@refallele.count/x@total.count
  }
  #random allele
  x.svd=ifelse(x.svd>matrix(runif(x@nsnp*npop),x@nsnp,npop),1,0)
 # ff[is.na(ff)]=NA #if freq=NA; test return NA also and scaling is done with na.rm=T
 # scaling (with dimension change to allow column operations)
  x.svd=scale(t(x.svd),scale=scale)
  x.svd[is.na(x.svd)]=0
  #SVD
  if(return.snploadings){nv=min(npop-1,x@nsnp-1)}else{nv=0}
  nu=min(npop-1,x@nsnp-1)
  x.svd=svd(x.svd,nu=nu,nv=nv)

  eigenvalues=x.svd$d[1:nu]**2 #eigenvalues (the last is 0 since centering)
  perc.var=100*eigenvalues/sum(eigenvalues)
  
  if(!is.null(plot.pcs)){
    i=plot.pcs[1] ; j=plot.pcs[2]
    plot(x.svd$u[,plot.pcs],
         xlab=paste0("PC",i," (",round(perc.var[i],2),"%)"),
         ylab=paste0("PC",j," (",round(perc.var[j],2),"%)"),...)
    text(x.svd$u[,plot.pcs],pos=3,popnames)
    abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
  }
  rownames(x.svd$u)=popnames
  if(return.snploadings){
    return(list(eigenvalues=eigenvalues,perc.var=perc.var,pop.loadings=x.svd$u,snp.loadings=x.svd$v))
  }else{
    return(list(eigenvalues=eigenvalues,perc.var=perc.var,pop.loadings=x.svd$u))
  }
  
}