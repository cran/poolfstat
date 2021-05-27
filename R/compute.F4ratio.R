#' Compute F4ratio (estimation of admixture rate) from an fstats object
#' @param x A fstats object containing estimates of fstats
#' @param num.quadruplet A character vector for the F4 quadruplet used in the F4ratio numerator (should be of the form "A,O;C,X" where A, O, C and X are the names of the population as defined in the countdata or pooldata object used to obtain fstats, see details)
#' @param den.quadruplet A character vector for the F4 quadruplet used in the F4ratio denominator (should be of the form "A,O;C,B" where A, O, C and B are the names of the populations as defined in the countdata or pooldata object used to obtain fstats, see details))
#' @details Assuming a 4 population phylogeny rooted with an outgroup O of the form (((A,B);C);O) and an admixed population X with two source populations related to B and C, the admixture rate alpha of the B-related ancestry is obtained using the ratio F4(A,O;C,X)/F4(A,O;C,B) (see Patterson et al., 2012 for more details).  
#' @return Either a scalar corresponding to the estimated admixture rate or, if F4 block-jackknife samples are available in the input fstats object (i.e., compute.fstats was run with return.F4.blockjackknife.samples = TRUE) a vector with three elements corresponding to the estimate of the admixture rate, the block-jacknife mean (may be slightly different than the previous since not exactly the same set of markers are used) and the standard error of the estimates.
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fstats=compute.fstats(pooldata)
#' @export

compute.f4ratio<-function(x,num.quadruplet,den.quadruplet){

if(!(is.fstats(x))){stop("Input data should either be an object of class fstats (see the function compute.fstats)\n")}
if(nrow(x@f4.values)==0){stop("fstats object must contain estimates of F4 (see compute.fstats that should be run with computeF4=TRUE)\n")}

f4.pop.names=x@comparisons[["F4"]] 

#recherhce des quadruplets
tmp.q=unlist(strsplit(num.quadruplet,split="[,;]"))
num.idx.f4=which(apply(f4.pop.names,1,f<-function(x){sum(x %in% tmp.q)==4 & ((sum(x[1:2] %in% tmp.q[1:2])==2) |(sum(x[1:2] %in% tmp.q[3:4])==2))}))
if(is.null(num.idx.f4)){
  stop("The num.quadruplet configuration is not available in the fstat object (check population names or fstats object)\n")
}

tmp.q=unlist(strsplit(den.quadruplet,split="[,;]"))
den.idx.f4=which(apply(f4.pop.names,1,f<-function(x){sum(x %in% tmp.q)==4 & ((sum(x[1:2] %in% tmp.q[1:2])==2) |(sum(x[1:2] %in% tmp.q[3:4])==2))}))
if(is.null(den.idx.f4)){
  stop("The den.quadruplet configuration is not available in the fstat object (check population names or fstats object)\n")
}


##check multiplicative factor
fact=1
available.f4=x@comparisons[["F4"]][num.idx.f4,] 
pair1=paste(available.f4[1:2],collapse=",") ; pair2=paste(available.f4[3:4],collapse=",")
tmp.q=unlist(strsplit(num.quadruplet,split="[,;]"))
target.pair1=paste(tmp.q[1:2],collapse=",") ; target.pair2=paste(tmp.q[3:4],collapse=",")
# if(pair1==target.pair1 & pair2!=target.pair2){fact=fact*-1}
# if(pair1!=target.pair1 & pair2==target.pair2){fact=fact*-1}
# if(pair2==target.pair1 & pair1!=target.pair2){fact=fact*-1}
# if(pair2!=target.pair1 & pair1==target.pair2){fact=fact*-1}
tst=(pair1==target.pair1) + (pair1==target.pair2) + (pair2==target.pair1) + (pair2==target.pair2)
if(tst==1){fact=-1*fact}

available.f4=x@comparisons[["F4"]][den.idx.f4,] 
pair1=paste(available.f4[1:2],collapse=",") ; pair2=paste(available.f4[3:4],collapse=",")
tmp.q=unlist(strsplit(den.quadruplet,split="[,;]"))
target.pair1=paste(tmp.q[1:2],collapse=",") ; target.pair2=paste(tmp.q[3:4],collapse=",")
tst=(pair1==target.pair1) + (pair1==target.pair2) + (pair2==target.pair1) + (pair2==target.pair2)
if(tst==1){fact=-1*fact}

#compute alpha

alpha=fact*x@f4.values[num.idx.f4,1]/x@f4.values[den.idx.f4,1]
n.blocks=ncol(x@F4.bjack.samples)
if(n.blocks>0){
  blocks.estimates=fact*x@F4.bjack.samples[num.idx.f4,]/x@F4.bjack.samples[den.idx.f4,]
  alpha.mean=mean(blocks.estimates)
  alpha.sd=sd(blocks.estimates)*(n.blocks-1)/sqrt(n.blocks)
  alpha=c(alpha,alpha.mean,alpha.sd)
  names(alpha)=c("Estimate","bjack mean","bjack s.e.")
}

return(alpha)
}


