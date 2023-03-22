#' Convert BayPass read count and haploid pool size input files into a pooldata object
#' @param genobaypass.file The name (or a path) of the BayPass read count file (see the BayPass manual \url{http://www1.montpellier.inra.fr/CBGP/software/baypass/})
#' @param poolsize.file The name (or a path) of the BayPass (haploid) pool size file (see the BayPass manual \url{http://www1.montpellier.inra.fr/CBGP/software/baypass/})
#' @param poolnames A character vector with the names of pool
#' @param snp.pos An optional two column matrix with nsnps rows containing the chromosome (or contig/scaffold) of origin and the position of each markers
#' @param min.cov.per.pool Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded
#' @param max.cov.per.pool Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overall read counts for the reference allele over the read coverage)
#' @param verbose If TRUE extra information is printed on the terminal
#' @details  Information on SNP position is only required for some graphical display or to carried out block-jacknife sampling estimation of confidence intervals. If no mapping information is given (default), SNPs will be assumed to be ordered on the same chromosome and separated by 1 bp. As blocks are defined with a number of consecutive SNPs (rather than a length), the latter assumption has actually no effect (except in the reported estimated block sizes in Mb).
#' @return A pooldata object containing 7 elements:
#' \enumerate{
#' \item "refallele.readcount": a matrix with nsnp rows and npools columns containing read counts for the reference allele (chosen arbitrarily) in each pool
#' \item "readcoverage": a matrix with nsnp rows and npools columns containing read coverage in each pool
#' \item "snp.info": a matrix with nsnp rows and four columns containing respectively the contig (or chromosome) name (1st column) and position (2nd column) of the SNP; the allele taken as reference in the refallele.readcount matrix (3rd column); and the alternative allele (4th column)
#' \item "poolsizes": a vector of length npools containing the haploid pool sizes
#' \item "poolnames": a vector of length npools containing the names of the pools
#' \item "nsnp": a scalar corresponding to the number of SNPs
#' \item "npools": a scalar corresponding to the number of pools
#' }
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  pooldata2genobaypass(pooldata=pooldata,writing.dir=tempdir())
#'  pooldata=genobaypass2pooldata(genobaypass.file=paste0(tempdir(),"/genobaypass"),
#'                                poolsize.file=paste0(tempdir(),"/poolsize"))
#' @export
genobaypass2pooldata<-function(genobaypass.file="",poolsize.file="",snp.pos=NA,poolnames=NA,min.cov.per.pool=-1,max.cov.per.pool=1e6,min.maf=-1,verbose=TRUE){
if(nchar(genobaypass.file)==0){stop("ERROR: Please provide the name of the read count data file (in BayPass format)")}
if(nchar(poolsize.file)==0){stop("ERROR: Please provide the name of the file containing haploid pool sizes (in BayPass format)")}
poolsizes=as.numeric(read.table(file(poolsize.file))[1,])

file.con=file(genobaypass.file,open="r") 
##
tmp.data=scan(file=file.con,nlines = 1,what="character",quiet=TRUE)
npools=length(tmp.data)/2
if(npools!=length(poolsizes)){stop("The number of pools in the Pool sizes file (number of elements of the first line) is not the same as the one in the Read count file (half the number of columns). Check both input files are in a valid BayPass format")}
close(file.con)
if(sum(is.na(poolnames))>0){
  poolnames=paste0("P",1:npools)
}else{
  poolnames=as.character(poolnames)
  if(length(poolnames)!=npools){stop("ERROR: The number of pools derived form the BayPass input files is different from the length of vector of pool names")}
}
###
pos.all2=(1:npools)*2 ; pos.all1=pos.all2-1
tmp.data=as.matrix(fread(file=genobaypass.file,data.table=FALSE))
#tmp.data=as.matrix(read.table(genobaypass.file,stringsAsFactors = F))
tmp.Y=tmp.data[,pos.all1]
tmp.N=tmp.data[,pos.all2]+tmp.Y
rm(tmp.data)
nsnps.tot=nrow(tmp.Y)
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
  poolnames=as.character(poolnames)
  if(length(poolnames)!=npools){stop("ERROR: The number of pools derived form the BayPass input files is different from the length of vector of pool names")}

tmp.maf=0.5-abs(0.5-rowSums(tmp.Y)/rowSums(tmp.N))
dum.sel=(rowSums(tmp.N>=min.cov.per.pool)==npools) & (rowSums(tmp.N<=max.cov.per.pool)==npools) & (tmp.maf>min.maf)

res<-new("pooldata")
res@npools=npools
res@nsnp=sum(dum.sel)
res@refallele.readcount=tmp.Y[dum.sel,]
rm(tmp.Y)
res@readcoverage=tmp.N[dum.sel,] 
rm(tmp.N)
tmp.snp.info=cbind(snp.pos[dum.sel,],matrix(NA,res@nsnp,2)) 
res@snp.info=data.frame(Chromosome=as.character(tmp.snp.info[,1]),
                        Position=as.numeric(tmp.snp.info[,2]),
                        RefAllele=as.character(tmp.snp.info[,3]),
                        AltAllele=as.character(tmp.snp.info[,4]),
                        stringsAsFactors = FALSE)
res@poolsizes=poolsizes
res@poolnames=poolnames

if(verbose){cat("Data consists of",res@nsnp,"SNPs for",res@npools,"Pools\n")}
return(res)
}
