#' Convert allele count input files from the Treemix program into a coundata object
#' @param genotreemix.file The name (or a path) of the Treemix allele count file (see the Treemix manual \url{https://bitbucket.org/nygcresearch/treemix/wiki/Home})
#' @param snp.pos An optional two column matrix with nsnps rows containing the chromosome (or contig/scaffold) of origin and the position of each markers
#' @param min.indgeno.per.pop  Minimal number of overall counts required in each population. If at least one pop is not genotyped for at least min.indgeno.per.pop (haploid) individual, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overall counts for the reference allele over the overall number of (haploid) individual genotyped)
#' @param verbose If TRUE extra information is printed on the terminal
#' @details  Information on SNP position is only required for some graphical display or to carried out block-jacknife sampling estimation of confidence intervals. If no mapping information is given (default), SNPs will be assumed to be ordered on the same chromosome and separated by 1 bp. As blocks are defined with a number of consecutive SNPs (rather than a length), the latter assumption has actually no effect (except in the reported estimated block sizes in Mb).
#' @return A countdata object containing 6 elements:
#' \enumerate{
#' \item "refallele.count": a matrix (nsnp rows and npops columns) with the allele counts for the reference allele
#' \item "total.count": a matrix (nsnp rows and npops columns) with the total number of counts (i.e., twice the number of genotyped individual for diploid species and autosomal markers)
#' \item "snp.info": a matrix with nsnp rows and four columns containing respectively the contig (or chromosome) name (1st column) and position (2nd column) of the SNP; the allele taken as reference in the refallele.count matrix (3rd column); and the alternative allele (4th column)
#' \item "popnames": a vector of length npops containing the names of the pops
#' \item "nsnp": a scalar corresponding to the number of SNPs
#' \item "npops": a scalar corresponding to the number of populations
#' }
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  ##NOTE: This example is just for the sake of illustration as it amounts 
#'  ##to interpret read count as allele count which must not be done in practice!
#'  dum=matrix(paste(pooldata@refallele.readcount,
#'    pooldata@readcoverage-pooldata@refallele.readcount,sep=","),
#'    ncol=pooldata@npools)
#'  colnames(dum)=pooldata@poolnames
#'  write.table(dum,file=paste0(tempdir(),"/genotreemix"),quote=FALSE,row.names=FALSE)
#'  countdata=genotreemix2countdata(genotreemix.file=paste0(tempdir(),"/genotreemix"))
#' @export
genotreemix2countdata<-function(genotreemix.file="",snp.pos=NA,min.indgeno.per.pop=-1,min.maf=-1,verbose=TRUE){
  if(nchar(genotreemix.file)==0){stop("ERROR: Please provide the name of the read count data file (in TreeMix format)")}
  if(verbose){cat("Begin reading data\n")}
  tmp.data=read.table(genotreemix.file,header=T,stringsAsFactors = F)  
  npops=ncol(tmp.data) ;  nsnps.tot=nrow(tmp.data)
  if(verbose){cat(npops,"populations and",nsnps.tot,"SNPs identified (parsing started)\n")}
  popnames=colnames(tmp.data)
  ###
  tmp.data=apply(as.matrix(tmp.data),2,f<-function(x){as.numeric(unlist(strsplit(x,split=",")))})
  pos.all2=(1:nsnps.tot)*2 ; pos.all1=pos.all2-1
  tmp.Y=tmp.data[pos.all1,]
  tmp.N=tmp.data[pos.all2,]+tmp.Y
  rm(tmp.data)
  if(sum(is.na(snp.pos))>0){
    snp.pos=cbind(rep("chr",nsnps.tot),1:nsnps.tot)
  }else{
    if(ncol(snp.pos)!=2 | nrow(snp.pos)!=nsnps.tot){
      cat("Warning: Invalid dimensions in snp.pos matrix\n")
      cat(" 2 columns expected and",ncol(snp.pos),"observed\n")
      cat(nsnps.tot," row expected (after genobaypass file) and",nrow(snp.pos),"observed\n") 
      cat("Default positions will be assumed (see details)\n")
      snp.pos=cbind(rep("chr",nsnps.tot),1:nsnps.tot)
    }else{
      snp.pos=as.matrix(snp.pos)
    }
  }
  tmp.maf=0.5-abs(0.5-rowSums(tmp.Y)/rowSums(tmp.N))
  dum.sel=(rowSums(tmp.N>=min.indgeno.per.pop)==npops) & (tmp.maf>min.maf)
  
  res<-new("countdata")
  res@npops=npops
  res@nsnp=sum(dum.sel)
  res@refallele.count=tmp.Y[dum.sel,]
  rm(tmp.Y)
  res@total.count=tmp.N[dum.sel,] 
  rm(tmp.N)
  tmp.snp.info=cbind(snp.pos[dum.sel,],matrix(NA,res@nsnp,2)) 
  res@snp.info=data.frame(Chromosome=as.character(tmp.snp.info[,1]),
                          Position=as.numeric(tmp.snp.info[,2]),
                          RefAllele=as.character(tmp.snp.info[,3]),
                          AltAllele=as.character(tmp.snp.info[,4]),
                          stringsAsFactors = FALSE)
  res@popnames=popnames
  colnames(res@refallele.count)=colnames(res@total.count)=res@popnames
  
  if(verbose){cat("The data set consists of",res@nsnp,"SNPs for",res@npops,"Pops\n")}
  return(res)
}
