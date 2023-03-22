#' Convert BayPass allele count input files into a coundata object
#' @param genobaypass.file The name (or a path) of the BayPass allele count file (see the BayPass manual \url{http://www1.montpellier.inra.fr/CBGP/software/baypass/})
#' @param popnames A character vector with the names of pool
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
#'  pooldata2genobaypass(pooldata=pooldata,writing.dir=tempdir())
#'  ##NOTE: This example is just for the sake of illustration as it amounts 
#'  ##to interpret read count as allele count which must not be done in practice!
#'  countdata=genobaypass2countdata(genobaypass.file=paste0(tempdir(),"/genobaypass")) 
#' @export
genobaypass2countdata<-function(genobaypass.file="",snp.pos=NA,popnames=NA,min.indgeno.per.pop=-1,min.maf=-1,verbose=TRUE){
  if(nchar(genobaypass.file)==0){stop("ERROR: Please provide the name of the read count data file (in BayPass format)")}
  tmp.data=as.matrix(fread(file=genobaypass.file,data.table=FALSE))
#  tmp.data=as.matrix(read.table(genobaypass.file,stringsAsFactors = F))  
  npops=ncol(tmp.data)/2
  nsnps.tot=nrow(tmp.data)
  cat(npops,"populations and",nsnps.tot,"SNPs read\n")
  if(sum(is.na(popnames))>0){
    popnames=paste0("P",1:npops)
  }else{
    popnames=as.character(popnames)
    if(length(popnames)!=npops){
      cat("Warning: the popnames vector length is",length(popnames),"while there is",npops,"in the BayPass genofile\n")
      cat("Default Name will be given to the populations\n")
      popnames=paste0("P",1:npops)
      }
  }
  ###
  pos.all2=(1:npops)*2 ; pos.all1=pos.all2-1
  tmp.Y=tmp.data[,pos.all1]
  tmp.N=tmp.data[,pos.all2]+tmp.Y
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
  
  if(verbose){cat("Data finally consists of",res@nsnp,"SNPs for",res@npops,"Pops\n")}
  return(res)
}
