#' Convert SelEstim read count input files into a pooldata object
#' @param genoselestim.file The name (or a path) of the SelEstim read count file (see the SelEstim manual \url{http://www1.montpellier.inra.fr/CBGP/software/selestim/})
#' @param poolnames A character vector with the names of pool
#' @param min.cov.per.pool Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded
#' @param max.cov.per.pool Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overal read counts for the reference allele over the read coverage)
#' @param nlines.per.readblock Number of Lines read simultaneously. Should be adapted to the available RAM.
#' @param verbose If TRUE extra information is printed on the terminal
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
#'  pooldata2genoselestim(pooldata=pooldata,writing.dir=tempdir())
#'  pooldata=genoselestim2pooldata(genoselestim.file=paste0(tempdir(),"/genoselestim"))
#' @export
genoselestim2pooldata<-function(genoselestim.file="",poolnames=NA,min.cov.per.pool=-1,max.cov.per.pool=1e6,min.maf=-1,nlines.per.readblock=1000000,verbose=TRUE){
if(nchar(genoselestim.file)==0){stop("ERROR: Please provide the name of the read count data file (in BayPass format)")}
file.con=file(genoselestim.file,open="r") 
##
tmp.data=scan(file=file.con,nlines = 2,what="character",quiet=TRUE)
npools=as.numeric(tmp.data[1])
tmp.data=scan(file=file.con,nlines = 1,what="character",quiet=TRUE)
poolsizes=as.numeric(tmp.data)
if(npools!=length(poolsizes)){stop("ERROR: problem in the number of pools or haploid pool size declaration. Check your selestim input file")}
if(sum(is.na(poolnames))>0){
  poolnames=paste0("Pool",1:npools)
}else{
  poolnames=as.character(poolnames)
  if(length(poolnames)!=npools){stop("ERROR: The number of pools derived form the BayPass input files is different from the length of vector of pool names")}
}
###
continue.reading=TRUE
nlines.read=0
time1=proc.time()
pos.all2=(1:npools)*2 ; pos.all1=pos.all2-1
while(continue.reading){
 tmp.data=matrix(as.numeric(scan(file=file.con,nlines = nlines.per.readblock,what="character",quiet=TRUE)),ncol=2*npools,byrow=T)  
 if(length(tmp.data)<nlines.per.readblock){continue.reading=FALSE}
 npos=nrow(tmp.data)
 if(npos>1){
  tmp.Y=tmp.data[,pos.all1] ; tmp.N=tmp.Y+tmp.data[,pos.all2]
  rm(tmp.data)
  ##filtres sur couverture et maf
  tmp.maf=0.5-abs(0.5-rowSums(tmp.Y)/rowSums(tmp.N))
  dum.sel=(rowSums(tmp.N>=min.cov.per.pool)==npools) & (rowSums(tmp.N<=max.cov.per.pool)==npools) & (tmp.maf>min.maf)
  tmp.Y=tmp.Y[dum.sel,] ; tmp.N=tmp.N[dum.sel,] 
  if(nlines.read==0){
   data.Y=tmp.Y ; data.N=tmp.N
  }else{
   data.Y=rbind(data.Y,tmp.Y) ; data.N=rbind(data.N,tmp.N)
 }
 nlines.read=nlines.read+npos
 if(verbose){
 cat(nlines.read/1000000,"millions lines processed in",round((proc.time()-time1)[1]/60,2)," min.; ",nrow(data.Y),"SNPs found\n")
 }
}
}
close(file.con)

res<-new("pooldata")
res@npools=npools
res@nsnp=nrow(data.Y)
res@refallele.readcount=data.Y
rm(data.Y)
res@readcoverage=data.N
rm(data.N)
tmp.snp.info=matrix(NA,res@nsnp,4) #q on rajoutera sliding windows, prevoir de changer en incluant la possibilite de fournir les positions
res@snp.info=data.frame(Chromosome=as.character(tmp.snp.info[,1]),
                        Position=as.numeric(tmp.snp.info[,2]),
                        RefAllele=as.character(tmp.snp.info[,3]),
                        AltAllele=as.character(tmp.snp.info[,4]),
                        stringsAsFactors = FALSE)
#rm(snpdet)
res@poolsizes=poolsizes
res@poolnames=poolnames

if(verbose){cat("Data consists of",res@nsnp,"SNPs for",res@npools,"Pools\n")}
return(res)
}
