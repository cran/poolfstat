#' Estimate the F-statistics (F2, F3, F3star, F4, Dstat) and within and across population diversity
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param nsnp.per.bjack.block Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling) 
#' @param computeDstat If TRUE compute Dstatistics (i.e. scaled F4). This may add some non negligible computation time if the number of population is large (n>15)
#' @param computeF3 If TRUE (default) compute all F3 and all F3star (i.e. scaled F3).
#' @param computeF4 If TRUE (default) compute all F4.
#' @param output.pairwise.fst If TRUE (default), output the npopxnpop matrix of pairwise-population Fst estimates (corresponding to the "Identity" method implemented in \code{\link{compute.pairwiseFST}}) in the pairwise.fst slot of the fstats output object (see help(fstats) for details) that may be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).
#' @param output.pairwise.div If TRUE (default), output the npopxnpop matrix of pairwise-population divergence (1-Q2) estimates  in the pairwise.div slot of the fstats output object (see help(fstats) for details) that may be visualized with e.g. heatmap function or used with a clustering function (e.g., hclust).
#' @param computeQmat If TRUE, compute the error covariance matrix between all F3 and F2 statistics (needed for admixture graph construction). This matrix may be very large if the number of pops is large. It is recommended to estimate it on a reduced sample of pops.
#' @param return.F2.blockjackknife.samples If TRUE (and nsnp.per.bjack.block>0) return an array of dimension (npopxnpopxnblocks) in an admixtools2 compatible format
#' @param return.F4.blockjackknife.samples Deprecated options (since v. 2.2.0) 
#' @param verbose If TRUE extra information is printed on the terminal
#' @details The function estimates for the n populations (or pools) represented in the input object x:
#' \enumerate{
#' \item The F2 statistics for all the \eqn{n(n-1)/2} pairs of populations (or pools) and their scaled version (equivalent, but faster, than Fst estimated with \code{\link{compute.pairwiseFST}} when method="Identity")
#' \item If n>2, The F3 statistics for all the \eqn{npools(npools-1)(npools-2)/2} possible triplets of populations (or pools) and their scaled version (named F3star after Patterson et al., 2012)
#' \item If n>3, The F4 statistics and the D-statistics (a scaled version of the F4) for all the \eqn{npools(npools-1)(npools-2)*(npools-3)/8} possible quadruplets of populations
#' \item The estimated within population heterozygosities (=1-Q1)
#' \item The estimated divergence for each pair of populations (=1-Q2) 
#' }
#' @return An object of class fstats (see help(fstats) for details)
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fstats=compute.fstats(pooldata)
#' @export
compute.fstats<-function(x,nsnp.per.bjack.block=0,computeDstat=FALSE,computeF3=TRUE,computeF4=TRUE,
                         output.pairwise.fst=TRUE,output.pairwise.div=TRUE,computeQmat=TRUE,
                         return.F2.blockjackknife.samples=FALSE,return.F4.blockjackknife.samples=FALSE,
                         verbose=TRUE){
  time1=proc.time()
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")
  } 
  if(is.pooldata(x)){npops=x@npools ; popnames=x@poolnames}
  if(is.countdata(x)){npops=x@npops ; popnames=x@popnames}  
  if(npops<4 & computeF4){computeF4=FALSE}#else{computeF4=TRUE}
  if(npops<3 & computeF3){computeF3=FALSE}#else{computeF3=TRUE}
  if(!computeF4){computeDstat=FALSE}
  if(!computeF3){computeQmat=FALSE}  
  if(return.F4.blockjackknife.samples){
    cat("NOTE: option return.F4.blockjackknife.samples has been inactivated since v. 2.2.0 of the package\n")
    cat("      option return.F2.blockjackknife.samples has been activated instead to allow computation of F4ratio\n")
    return.F2.blockjackknife.samples=TRUE
  }
  
  out=new("fstats")
  out@comparisons[["pops"]]=popnames
  #########################
  #####jacknife: if not then res.size=1 else res.size=nbloc+1 (this is the dimension of the Q1 et Q2 per bloc estimates corresponding to the sum of within bloc SNP-specific sum; the last dimension is used for the overall sum of SNP-specific stat)
  if(nsnp.per.bjack.block>0){out@blockjacknife=TRUE}else{out@blockjacknife=FALSE}
  if(!out@blockjacknife){
    if(return.F2.blockjackknife.samples){stop("Incompatible Option: nsnp.per.bjack.block must be >0 if return.F2.blockjackknife.samples is TRUE\n")}   
    computeQmat=FALSE #stop("Incompatible Option: nsnp.per.bjack.block must be >0 if computeQmat is TRUE\n")  
    block.id=rep(0,x@nsnp)
    #    block.id[is.na(bjack.blocks$snp.block.id)]=0
    nblocks=0   
  }else{
    if(verbose){cat("Block-Jackknife specification\n")}
    bjack.blocks=generate.jackknife.blocks(x,nsnp.per.bjack.block,verbose=verbose)
    block.id=bjack.blocks$snp.block.id
    block.id[is.na(bjack.blocks$snp.block.id)]=0
    nblocks=bjack.blocks$nblocks
    rm(bjack.blocks)
  }
  #########################  
  #########################  
  
  #Computing Q1 for each pop  
  if(verbose){cat("Estimating Q1\n")}
  if(is.pooldata(x)){#verbose inactive dans calul cpp car si peu de pops ca peut bugger
    Q1=.compute_H1(refcount = x@refallele.readcount,totcount = x@readcoverage,nblocks = nblocks,block_id = block.id-1,verbose=verbose)
    Q1=1-Q1*x@poolsizes/(x@poolsizes-1)
  }else{
    Q1=1-.compute_H1(refcount = x@refallele.count,totcount = x@total.count,nblocks = nblocks,block_id = block.id-1,verbose=verbose)
  }
  #Computing Q2 for each snp/pop pair
  if(verbose){cat("Estimating Q2\n")}
  if(is.pooldata(x)){
    Q2=.compute_Q2(refcount = x@refallele.readcount,totcount = x@readcoverage,nblocks = nblocks,block_id = block.id-1,verbose=verbose)
  }else{
    Q2=.compute_Q2(refcount = x@refallele.count,totcount = x@total.count,nblocks = nblocks,block_id = block.id-1,verbose=verbose) 
  }
  #Computing F2
  n.f2=npops*(npops-1)/2 #nrow(Q2)
  p1p2.idx=combn(npops,2)   #index=((tmp.cb[1,]-1)*npops + tmp.cb[2,])-(tmp.cb[1,]*(tmp.cb[1,]+1)/2)
  tmp.names=cbind(popnames[p1p2.idx[1,]],popnames[p1p2.idx[2,]])
  colnames(tmp.names)=c("P1","P2")
  rownames(tmp.names)=paste(tmp.names[,1],tmp.names[,2],sep=",")
  out@comparisons[["F2"]]=out@comparisons[["Fst"]]=tmp.names
  rm(tmp.names)
  F2=(Q1[p1p2.idx[1,],]+Q1[p1p2.idx[2,],])/2 - Q2
  ############
  ##if blockjackknige compute l.o.o. value for Q1, Q2 and F2 (and store them in the original objects) to derive other stats
  ############
  if(out@blockjacknife){
    if(return.F2.blockjackknife.samples){
      #export F2 in admixtools2 format before computing l.o.o.
      out@F2.bjack.samples=array(0,dim=c(npops,npops,nblocks),
                                 dimnames=list(popnames,popnames,rep(paste0("l",nsnp.per.bjack.block),nblocks)))
      dum.mat=matrix(0,npops,npops) ; dum.up=upper.tri(dum.mat) ; dum.lo=lower.tri(dum.mat)
      for(i in 1:nblocks){
        dum.mat[dum.lo]=F2[,i] ; dum.mat[dum.up]=t(dum.mat)[dum.up] ; out@F2.bjack.samples[,,i]=dum.mat
      }
      rm(dum.mat)
    }
    #
    Q1[,1:nblocks]=(rowSums(Q1[,1:nblocks]) - Q1[,1:nblocks])/(nblocks - 1)
    Q2[,1:nblocks]=(rowSums(Q2[,1:nblocks]) - Q2[,1:nblocks])/(nblocks - 1)   
    F2[,1:nblocks]=(rowSums(F2[,1:nblocks]) - F2[,1:nblocks])/(nblocks - 1)  
    sd.corr=(nblocks-1)/sqrt(nblocks)
  }
  
  ###Computing heterozygosity
  if(verbose){cat("Estimating within-population heterozygosities\n")}
  if(out@blockjacknife){
    out@heterozygosities=as.data.frame(cbind(1-Q1[,nblocks+1],rowMeans(1-Q1[,1:nblocks]),apply(1-Q1[,1:nblocks],1,sd)*sd.corr))
    colnames(out@heterozygosities)=c("Estimate","bjack mean","bjack s.e.")
    rownames(out@heterozygosities)=popnames
  }else{
    out@heterozygosities=as.data.frame(matrix(1-Q1,ncol=1,dimnames = list(popnames,"Estimate")))
  }
  ###Computing F2
  if(verbose){cat("Estimating F2, pairwise Fst (method=Identity), and pairwise divergence (1-Q2)\n")}
  if(out@blockjacknife){
    out@divergence=as.data.frame(cbind(1-Q2[,nblocks+1],rowMeans(1-Q2[,1:nblocks]),apply(1-Q2[,1:nblocks],1,sd)*sd.corr))
    out@f2.values=as.data.frame(cbind(F2[,nblocks+1],rowMeans(F2[,1:nblocks]),apply(F2[,1:nblocks],1,sd)*sd.corr))
    tmp=F2/(1-Q2)
    out@fst.values=as.data.frame(cbind(tmp[,nblocks+1],rowMeans(tmp[,1:nblocks]),apply(tmp[,1:nblocks],1,sd)*sd.corr))   
    colnames(out@f2.values)=colnames(out@fst.values)=colnames(out@divergence)=c("Estimate","bjack mean","bjack s.e.")
    rownames(out@f2.values)=rownames(out@fst.values)=rownames(out@divergence)=rownames(out@comparisons[["F2"]])
    rm(tmp) ; gc()
  }else{
    out@divergence=as.data.frame(matrix(1-Q2,ncol=1,dimnames = list(rownames(out@comparisons[["F2"]]),"Estimate")))
    out@f2.values=as.data.frame(matrix(F2,ncol=1,dimnames = list(rownames(out@comparisons[["F2"]]),"Estimate")))
    out@fst.values=as.data.frame(matrix(F2/(1-Q2),ncol=1,dimnames = list(rownames(out@comparisons[["F2"]]),"Estimate"))) 
  }
  if(output.pairwise.fst){
    out@pairwise.fst=matrix(0,npops,npops,dimnames=list(popnames,popnames))
    out@pairwise.fst[lower.tri(out@pairwise.fst)]=out@fst.values$Estimate
    out@pairwise.fst[upper.tri(out@pairwise.fst)]=t(out@pairwise.fst)[upper.tri(out@pairwise.fst)]    
  }
  if(output.pairwise.div){
    out@pairwise.div=matrix(0,npops,npops,dimnames=list(popnames,popnames))
    out@pairwise.div[lower.tri(out@pairwise.div)]=out@divergence$Estimate
    out@pairwise.div[upper.tri(out@pairwise.div)]=t(out@pairwise.div)[upper.tri(out@pairwise.div)]    
  }
  ###Compute F3 from F2: F3(A;B,C)=(F2(A,B)+F2(A,C)-F2(B,C))/2
  if(computeF3){
    n.f3=npops*(npops-1)*(npops-2)/2 #npop*choose(npop-1,2)
    if(verbose){cat("Estimating F3 and F3* (n=",n.f3,"configurations)\n")}
    out@comparisons[["F3"]]=.generateF3names(popnames)
    rownames(out@comparisons[["F3"]])=out@comparisons[["F3"]][,1]
    out@comparisons[["F3"]]=out@comparisons[["F3"]][,-1]
    colnames(out@comparisons[["F3"]])=c("Px","P1","P2")
    out@comparisons[["F3star"]]=out@comparisons[["F3"]]
    if(out@blockjacknife){
      tmp.est=.compute_F3fromF2(F2val = F2[,nblocks+1],Hval = 1-Q1[,nblocks+1],npops = npops)
      tmp.bjack=.compute_F3fromF2samples(blockF2=F2[,1:nblocks],blockHet = 1-Q1[,1:nblocks],npops=npops,verbose = verbose)
      out@f3.values=as.data.frame(cbind(tmp.est[,1],tmp.bjack[,1],tmp.bjack[,2],tmp.bjack[,1]/tmp.bjack[,2]))
      out@f3star.values=as.data.frame(cbind(tmp.est[,2],tmp.bjack[,3],tmp.bjack[,4],tmp.bjack[,3]/tmp.bjack[,4]))      
      rownames(out@f3.values)=rownames(out@f3star.values)=rownames(out@comparisons[["F3"]])
      colnames(out@f3.values)=colnames(out@f3star.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")    
      rm(tmp.bjack)
    }else{
      tmp.est=.compute_F3fromF2(F2val = F2,Hval = 1-Q1,npops = npops)
      out@f3.values=as.data.frame(matrix(tmp.est[,1],ncol=1,dimnames = list(rownames(out@comparisons[["F3"]]),"Estimate")))
      out@f3star.values=as.data.frame(matrix(tmp.est[,2],ncol=1,dimnames = list(rownames(out@comparisons[["F3"]]),"Estimate")))      
    }
    rm(tmp.est) ; gc()
  }
  ###compute F4 and D from F2: F4(A,B;C,D)=(F2(A,D)+F2(B,C)-F2(A,C)-F2(B,D))/2
  ###cf function   .compute_F4fromF2(F2val=F2[,nblocks+1],npops = npops) si on veut shunter Dstat (est ce necessaire maintenant?)
  if(computeF4){
    n.f4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 #choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
    out@comparisons[["F4"]]=.generateF4names(popnames)
    rownames(out@comparisons[["F4"]])=out@comparisons[["F4"]][,1]
    out@comparisons[["F4"]]=out@comparisons[["F4"]][,-1]
    colnames(out@comparisons[["F4"]])=c("P1","P2","P3","P4")
    if(computeDstat){
      if(verbose){cat("Estimating F4 and Dstat (n=",n.f4,"configurations)\n")}
      if(verbose){cat("  Step 1/2: Estimating Dstat denominators\n")}      
      if(is.pooldata(x)){
        tmp.denom=.compute_blockDdenom(refcount = x@refallele.readcount,totcount = x@readcoverage,nblocks = nblocks,block_id = block.id-1,verbose=verbose)
      }else{
        tmp.denom=.compute_blockDdenom(refcount = x@refallele.count,totcount = x@total.count,nblocks = nblocks,block_id = block.id-1,verbose=verbose) 
      }
      if(out@blockjacknife){
        tmp.denom[,1:nblocks]=(rowSums(tmp.denom[,1:nblocks]) - tmp.denom[,1:nblocks])/(nblocks - 1) #l.o.o estimates for blocks
      }
      out@comparisons[["Dstat"]]=out@comparisons[["Dstat"]]
    }else{
      if(verbose){cat("Estimating F4 (n=",n.f4,"configurations)\n")}
    }
    
    tmp.est=.compute_F4fromF2(F2val = F2[,nblocks+1],npops = npops)
    if(out@blockjacknife){
      if(computeDstat){
        if(verbose){cat("  Step 2/2: Estimating F4 and D- statistics\n")}   
        tmp.bjack=.compute_F4DfromF2samples(blockF2=F2[,1:nblocks],blockDenom = tmp.denom,npops=npops,verbose = verbose)
        out@f4.values=as.data.frame(cbind(tmp.est,tmp.bjack[,1],tmp.bjack[,2],tmp.bjack[,1]/tmp.bjack[,2]))
        out@Dstat.values=as.data.frame(cbind(tmp.est/tmp.denom[,nblocks+1],tmp.bjack[,3],tmp.bjack[,4],tmp.bjack[,3]/tmp.bjack[,4]))
        rownames(out@f4.values)=rownames(out@Dstat.values)=rownames(out@comparisons[["F4"]])
        colnames(out@f4.values)=colnames(out@Dstat.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")   
      }else{
        tmp.bjack=.compute_F4fromF2samples(blockF2=F2[,1:nblocks],npops=npops,verbose = verbose)
        out@f4.values=as.data.frame(cbind(tmp.est,tmp.bjack[,1],tmp.bjack[,2],tmp.bjack[,1]/tmp.bjack[,2]))
        rownames(out@f4.values)=rownames(out@comparisons[["F4"]])
        colnames(out@f4.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")           
      }
      rm(tmp.bjack)
    }else{
      out@f4.values=as.data.frame(matrix(tmp.est,ncol=1,dimnames = list(rownames(out@comparisons[["F4"]]),"Estimate")))
      if(computeDstat){
        out@Dstat.values=as.data.frame(matrix(tmp.est/tmp.denom,ncol=1,dimnames = list(rownames(out@comparisons[["F4"]]),"Estimate")))     
      }
    }    
    rm(tmp.est) ; gc()
  }
  ################
  ################
  if(computeQmat){
    if(verbose){cat("Estimating Qmat, the error covariance matrix (",n.f2+n.f3,"x",n.f2+n.f3,")\n")}
    try(out@Q.matrix<-.compute_QmatfromF2samples(blockF2 = F2[,-(nblocks+1)],npops = npops,verbose=verbose),silent=TRUE)
    if(nrow(out@Q.matrix)==0){
      cat("Warning: not enough memory to compute and store the error covariance matrix (Q.matrix).\nIt will not be possible to fit graph with the output fstats object\nYou may try to run compute.fstats on a smaller object, i.e., with less populations\n(use pooldata.subset or countdata.subset to create it) involving only the populations of interest.\n")
    }else{
      rownames(out@Q.matrix)=colnames(out@Q.matrix)=c(rownames(out@comparisons[["F2"]]),rownames(out@comparisons[["F3"]]))
    }
  }
  
  if(verbose){
    time.elapsed=((proc.time()-time1)[3])
    nhours=floor(time.elapsed/3600) ; nminutes=floor((time.elapsed-nhours*3600)/60) ;  nseconds=round(time.elapsed-nhours*3600-nminutes*60)  
    cat("\nOverall Analysis Time:",nhours,"h",nminutes, "m",nseconds,"s\n")  
  }
  
  return(out)
  
}
