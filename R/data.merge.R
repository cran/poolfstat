#' Merge two pooldata or countdata objects
#' @param x1 First pooldata or countdata object to merge
#' @param x2 Second pooldata or countdata object to merge
#' @param fake.pool.size Specifies the haploid sample size used when merging a `countdata` object with a `pooldata` object to create a pseudo pooldata object containing all samples (default = 1e6), see details.
#' @param verbose If TRUE return some information
#' @details This function merges two objects of class `pooldata` and/or `countdata`, automatically checking their structure for consistency.
#' The merging behavior depends on the relationship between sample names and SNP identifiers:
#' 
#' \strong{1. Merging different samples (same SNPs):} If SNP names are identical but (pool or population) sample names differ, the function merges data from the distinct samples into a single `pooldata` or `countdata` object that includes all samples.
#' 
#' \strong{2. Merging different SNPs (same sample):} If sample names are identical but SNP names differ, the SNP data from each object are merged for each shared sample, effectively combining the variant information into one object.
#' 
#' \strong{3. Merging a `countdata` object with a `pooldata` object:} In this case, the function returns a `pooldata` object. Allele counts from the `countdata` object are converted into pseudo read counts. 
#'   To ensure compatibility, the haploid sample size for the sample originally contained in the `countdata` object is set to the value specified by the `fake.haploid.size` argument (default = \code{1e6}).
#'   Setting this value to a very large number (as in the default) ensures that each read count is treated as originating from a distinct haploid individual—
#'   mimicking Pool-Seq data where read coverage is much lower than the haploid sample size. This effectively disables Pool-Seq-specific bias corrections in downstream statistical analyses.
#'   Importantly, when merging objects of different types, only SNP-level merging is permitted. In this context, population samples are indeed expected to be necessarily distinct (at least in terms of effective haploid sample sizes).
#'   
#' @return A new `pooldata` or `countdata` object, depending on the input types.
#' @seealso To obtain description of the `countdata` and `pooldata` objects, see \code{\link{countdata}} and \code{\link{pooldata}}
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata1=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  pooldata2=pooldata1
#'  #Merge pooldata1 and pooldata2 by SNP
#'  pooldata2@poolnames=paste0(pooldata2@poolnames,"_2") #pool names must be different
#'  data.merged=data.merge(pooldata1,pooldata2)
#'  #Merge pooldata1 and pooldata2 by POP
#'  pooldata2=pooldata1
#'  pooldata2@snp.info[,1]=paste0(pooldata2@snp.info[,1],"_2") #SNP info must be different
#'  data.merged=data.merge(pooldata1,pooldata2)  
#'  #Merge pooldata1 with a countdata object
#'  #create a countdata object (NOTE: This example is just for the sake of illustration)
#'  pooldata2genobaypass(pooldata=pooldata1,writing.dir=tempdir())
#'  countdata=genobaypass2countdata(genobaypass.file=paste0(tempdir(),"/genobaypass")) 
#'  countdata@snp.info=pooldata1@snp.info
#'  countdata@popnames=paste0(countdata@popnames,"_2") #pop names must be different
#'  data.merged=data.merge(pooldata1,countdata)  
#' @export

data.merge<-function(x1,x2,fake.pool.size=1e6,verbose=TRUE){
  if(!(is.pooldata(x1)) & !(is.countdata(x1))) {stop("Object x1 is not a valid pooldata or countdata object\n)")}
  if(!(is.pooldata(x2)) & !(is.countdata(x2))) {stop("Object x2 is not a valid pooldata or countdata object\n)")}
  if(is.pooldata(x1)){samp.name.1=x1@poolnames}else{samp.name.1=x1@popnames}
  if(is.pooldata(x2)){samp.name.2=x2@poolnames}else{samp.name.2=x2@popnames}  
  mergeSNP=mergePOP=FALSE
  if(identical(samp.name.1,samp.name.2)){mergeSNP=TRUE}
  if(identical(x1@snp.info[,1:2],x2@snp.info[,1:2])){mergePOP=TRUE} 
  if(mergeSNP+mergePOP!=1){stop("Objects x1 and x2 cannot be merged because i) SNP information (chromosome/position) or sample names are inconsistent (i.e., at least one common id); or ii) are both the same\n")}
  if(is.pooldata(x1) & is.pooldata(x2)){#both pooldata
    out=x1
    if(mergePOP){
      out@npools=out@npools+x2@npools
      out@refallele.readcount=cbind(out@refallele.readcount,x2@refallele.readcount)
      out@readcoverage=cbind(out@readcoverage,x2@readcoverage)    
      out@poolsizes=c(out@poolsizes,x2@poolsizes)
      out@poolnames=c(out@poolnames,x2@poolnames)      
    }
    if(mergeSNP){
      out@nsnp=out@nsnp+x2@nsnp
      out@refallele.readcount=rbind(out@refallele.readcount,x2@refallele.readcount)
      out@readcoverage=rbind(out@readcoverage,x2@readcoverage)    
      out@snp.info=rbind(out@snp.info,x2@snp.info)
    }
  }else{  
   if(is.countdata(x1) & is.countdata(x2)){#both countdata
     out=x1
     if(mergePOP){
      out@npops=out@npops+x2@npops
      out@refallele.count=cbind(out@refallele.count,x2@refallele.count)
      out@total.count=cbind(out@total.count,x2@total.count)    
      out@popnames=c(out@popnames,x2@popnames)      
     }
     if(mergeSNP){
       out@nsnp=out@nsnp+x2@nsnp
       out@refallele.count=rbind(out@refallele.count,x2@refallele.count)
       out@total.count=rbind(out@total.count,x2@total.count)    
       out@snp.info=rbind(out@snp.info,x2@snp.info)
     }
   }else{#countdata and pooldata (can only be merged by sample name)
     if(mergeSNP){
       stop("The objects x1 ans x2 cannot be merged by SNP because of different types. If you really want to do this, change pop names of the count data object to specify it is a different sample\n")
     }   
     if(is.countdata(x1)){
       out=x2
       out@npools=out@npools+x1@npops
       out@refallele.readcount=cbind(out@refallele.readcount,x1@refallele.count)
       out@readcoverage=cbind(out@readcoverage,x1@total.count)    
       out@poolsizes=c(out@poolsizes,rep(fake.pool.size,x1@npops))
       out@poolnames=c(out@poolnames,x1@popnames)     
     }else{
       out=x1
       out@npools=out@npools+x2@npops
       out@refallele.readcount=cbind(out@refallele.readcount,x2@refallele.count)
       out@readcoverage=cbind(out@readcoverage,x2@total.count)    
       out@poolsizes=c(out@poolsizes,rep(fake.pool.size,x2@npops))
       out@poolnames=c(out@poolnames,x2@popnames)           
     }
     }
  }
  if(verbose){
    cat("The two objects have been successfully merged.\n")
    if(is.countdata(out)){
     cat("The merged allele count data set consists of",out@nsnp,"SNPs for",out@npops,"Pop samples\n")
    }else{
      cat("The merged (PoolSeq) read count data set consists of",out@nsnp,"SNPs for",out@npools,"Pooled samples\n")      
    }  
  }
  return(out)
}
