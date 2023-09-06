#' Compute F4ratio (estimation of admixture rate) from an fstats object
#' @param x A fstats object containing estimates of fstats
#' @param num.quadruplet A character vector for the F4 quadruplet used in the F4ratio numerator (should be of the form "A,O;C,X" where A, O, C and X are the names of the population as defined in the countdata or pooldata object used to obtain fstats, see details)
#' @param den.quadruplet A character vector for the F4 quadruplet used in the F4ratio denominator (should be of the form "A,O;C,B" where A, O, C and B are the names of the populations as defined in the countdata or pooldata object used to obtain fstats, see details))
#' @details Assuming a 4 population phylogeny rooted with an outgroup O of the form (((A,B);C);O) and an admixed population X with two source populations related to B and C, the admixture rate alpha of the B-related ancestry is obtained using the ratio F4(A,O;C,X)/F4(A,O;C,B) (see Patterson et al., 2012 for more details).  
#' @return Either a scalar corresponding to the estimated admixture rate or, if F2 block-jackknife samples are available in the input fstats object (i.e., compute.fstats was run with return.F2.blockjackknife.samples = TRUE) a vector with three elements corresponding to the estimate of the admixture rate, the block-jacknife mean (may be slightly different than the previous since not exactly the same set of markers are used) and the standard error of the estimates.
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fstats=compute.fstats(pooldata)
#' @export

compute.f4ratio<-function(x,num.quadruplet,den.quadruplet){
  
  if(!(is.fstats(x))){stop("Input data should either be an object of class fstats (see the function compute.fstats)\n")}
  if(nrow(x@f2.values)==0){stop("fstats object must contain estimates of F4 (see compute.fstats that should be run with computeF4=TRUE)\n")}
  
  #find index of the different pops (to compute F4 from F2 and to find F2 in the array)
  num.pops=unlist(strsplit(num.quadruplet,split="[,;]"))  
  if(length(num.pops)!=4){
    stop("The num.quadruplet configuration could not be parsed (check pop. names or configuration name)\n")
  }
  den.pops=unlist(strsplit(den.quadruplet,split="[,;]"))    
  if(length(den.pops)!=4){
    stop("The den.quadruplet configuration could not be parsed (check pop. names or configuration name)\n")
  }
  
  num.pops.idx=den.pops.idx=rep(NA,4)
  for(i in 1:4){
    tmp=which(x@comparisons$pops==num.pops[i])
    if(length(tmp)>0){num.pops.idx[i]=tmp}else{
      stop("Pop. named ",num.pops[i]," from den.quadruplet configuration was not included in the fstats object analysis\n") 
    }
    tmp=which(x@comparisons$pops==den.pops[i])
    if(length(tmp)>0){den.pops.idx[i]=tmp}else{
      stop("Pop. named ",den.pops[i]," from den.quadruplet configuration was not included in the fstats object analysis\n")  
    }
  }
  
  #f2 idx for F2(A,D), F2(B,C), F2(A,C), and F2(B,D) required to compute F4
  f2.required=rbind(c(1,4),c(2,3),c(1,3),c(2,4))
  npops.tot=length(x@comparisons$pops) #total number of pop in the fstats object
  tmp=cbind(num.pops.idx[f2.required[,1]],num.pops.idx[f2.required[,2]])
  tmp=cbind(apply(tmp,1,min),apply(tmp,1,max))
  num.f2.idx=(tmp[,1]-1)*npops.tot + tmp[,2] - tmp[,1]*(tmp[,1]+1)/2
  tmp=cbind(den.pops.idx[f2.required[,1]],den.pops.idx[f2.required[,2]])
  tmp=cbind(apply(tmp,1,min),apply(tmp,1,max))
  den.f2.idx=(tmp[,1]-1)*npops.tot + tmp[,2] - tmp[,1]*(tmp[,1]+1)/2
  
  #compute alpha
  num.f4=x@f2.values$Estimate[num.f2.idx[1]] + x@f2.values$Estimate[num.f2.idx[2]] - 
    x@f2.values$Estimate[num.f2.idx[3]] - x@f2.values$Estimate[num.f2.idx[4]] 
  den.f4=x@f2.values$Estimate[den.f2.idx[1]] + x@f2.values$Estimate[den.f2.idx[2]] - 
    x@f2.values$Estimate[den.f2.idx[3]] - x@f2.values$Estimate[den.f2.idx[4]] 
  
  alpha=num.f4/den.f4
  
  if(length(x@F2.bjack.samples)>0){
    n.blocks=dim(x@F2.bjack.samples)[3]  
    num.f2.bjack=den.f2.bjack=matrix(0,4,n.blocks)
    for(i in 1:4){
      num.f2.bjack[i,]=x@F2.bjack.samples[num.pops.idx[f2.required[i,1]],num.pops.idx[f2.required[i,2]],]
      den.f2.bjack[i,]=x@F2.bjack.samples[den.pops.idx[f2.required[i,1]],den.pops.idx[f2.required[i,2]],]      
    }
    #compute l.o.o
    num.f2.bjack=(rowSums(num.f2.bjack) - num.f2.bjack)/(n.blocks - 1)  
    den.f2.bjack=(rowSums(den.f2.bjack) - den.f2.bjack)/(n.blocks - 1) 
    bj.f4rat=(colSums(num.f2.bjack[1:2,])-colSums(num.f2.bjack[3:4,]))/
      (colSums(den.f2.bjack[1:2,])-colSums(den.f2.bjack[3:4,]))
    alpha.mean=mean(bj.f4rat)
    alpha.se=sd(bj.f4rat)*(n.blocks-1)/sqrt(n.blocks)
    alpha=c(alpha,alpha.mean,alpha.se)
    names(alpha)=c("Estimate","bjack mean","bjack s.e.")
  }
  
  return(alpha)
}


