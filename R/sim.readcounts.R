#' Simulate read counts from count data and return a pooldata object
#' @param x A countdata object containing allele count information
#' @param lambda.cov Numeric vector of length npop giving the expected coverage of each pool
#' @param overdisp Numeric value giving overdispersion of coverages (see details)
#' @param min.rc Integer giving the minimal read count for an allele to be considered as true allele
#' @param maf.thr Float giving the MAF threshold for SNP filtering
#' @param seq.eps Numeric value giving the sequencing error rate
#' @param exp.eps Numeric value giving the experimental error leading to unequal contribution of individual to the pool reads
#' @param genome.size Size of the genome (only considered when seq.eps>0 to simulated spurious SNPs generated at monomorphic position)
#' @param verbose If TRUE extra information is printed on the terminal
#' @details The function implements a simulation approach similar to that described in Gautier et al. (2021). 
#' Read coverages are sampled from a distribution specified by the lambda.cov vector and the overdisp scalar. Note that overdisp is the same for all pop sample but the expected coveragese (specified in the lambda.cov vector) may vary across pool. If overdisp=1 (default), coverages are assumed Poisson distributed with mean (and variance) equal to the value specified in the lambda.cov vector. If overdisp>1, coverages follows a Negative Binomial distribution with a mean equal to lambda and variance equal to overdisp*lambda. Finally, if overdisp<1, no variation in coverage is introduced and all coverages are equal to the value specified in the lambda vector although they may (slightly) vary in the output when seq.eps>0 due to the removal of error reads.
#' The seq.eps parameter control sequencing error rate. Sequencing errors are modeled following Gautier et al. (2021) i.e. read counts for the four possible bases are sampled from a multinomial distribution Multinom(c,\{f*(1-eps)+(1-f)*eps/3;f*eps/3+(1-f)*(1-eps),eps/3,eps/3\})  where c is the read coverage and f the reference allele frequencies (obtained from the count data). When seq.eps>0, spurious SNPs may be generated at monomorphic positions (the number of which being equal to the size of the genome, provided with the genome.size argument, minus the number of SNPs in the countdata object). These spurious SNPs are simulated using the same error model (Multinom(c,\{1-eps,eps/3,eps/3,eps/3\}). Only bi-allelic SNPs passing filtering conditions specified by min.rc (which controls the minimal read count for an allele to be deemed as true, i.e. if more than two alleles have >= min.rc counts then the SNP is excluded because non-bi-allelic) and maf.thr (threshold on the major allele frequency computed over all read counts) are included in the output.
#'  Experimental error exp.eps control the contribution of individual (assumed diploid) to the pools following the model described in Gautier et al. (2013). The parameter exp.eps corresponds to the coefficient of variation of the individual contributions. For example, in a pool of 10 individuals and a Poisson distributed coverage of mean 100, exp.eps=0.5 correspond to a situation where the 5 most contributing individuals contribute $>2$ times reads than the others. When exp.eps tends toward 0, all individuals contribute equally to the pool and there is no experimental error. Note that the number of (diploid) individuals for each SNP and pop. sample is deduced from the input total count   (it may thus differ over SNP when the total counts are not the same). 
#' @return A pooldata object containing simulated read counts
#' @seealso To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  #not run 
#' @export
sim.readcounts<-function(x,lambda.cov=rep(50,x@npops),overdisp=1,seq.eps=0,exp.eps=0,maf.thr=0,min.rc=2,genome.size=0,verbose=TRUE){
  if(!is.countdata(x)){stop("x must be a countdata object\n")}
  npop=x@npops 
  if(length(lambda.cov)!=npop){stop("lambda.cov must be a vector of same length as teh number of population in the input countdata object\n")}
  if(seq.eps<0 | seq.eps>1){stop("seq.eps must be between 0 and 1\n")}
  if(exp.eps<0){stop("exp.eps must be >0\n")}  
  out=new("pooldata")
  out@poolnames=x@popnames ;  out@poolsizes=apply(x@total.count,2,max) ; out@npools=npop
  ##simu pos polymorphe
  if(verbose){cat("Start simulation of Read counts for the",x@nsnp,"polymorphic SNPs included in the count data object\n")}
  tmp=.simureads_poly(x@refallele.count,x@total.count,lambda.cov,overdisp,min.rc,maf.thr,seq.eps,exp.eps)
  pos.poly=rowSums(tmp)!=0
  if(sum(pos.poly)>0){
   out@refallele.readcount=tmp[pos.poly,1:npop]
   out@readcoverage=tmp[pos.poly,-1*(1:npop)]
   out@snp.info=x@snp.info[pos.poly,]
   out@snp.info[,3]="poly"
  }else{
    stop("No SNP left among the polymorophic ones: seq.eps may be too high (generating too many multi-allelic SNPs) or SNP filtering criteria (min.maf and min.rc) too stringent\n")
  }
   rm(tmp)
  if(verbose){cat("Simulation ended:",sum(pos.poly),"SNPs retained\n")}
  if(seq.eps>0){
    npos=floor(genome.size-sum(pos.poly))
    if(npos<1){cat("Warning: no simulation of spurious SNPs in monomorphic position because the genome size (specified genome.size argument) is smaller than the number of polymorphic SNPs\n")
      }else{
        if(verbose){cat("Start simulation of spurious SNPs for",npos,"monomorphic positions\n")}
        tmp=.simureads_mono(npos,npop,lambda.cov,overdisp,min.rc,maf.thr,seq.eps)
        nsnperr=nrow(tmp)
        if(verbose){cat("Simulation of spurious SNPs ended:",nsnperr,"SNPs retained\n")}
        if(nsnperr>0){
          chr.start=table(x@snp.info[,1]) #initilise to number of SNPs per chr
          chr.names=names(chr.start) ; nchr=length(chr.names)
          chr.end   = ceiling(cumsum(chr.start*genome.size/sum(chr.start))) #chr length in proportion of initial SNPs
          chr.start = c(1,chr.end[-nchr]-1)
          snperr.coord=sample(1:genome.size,nsnperr,replace=TRUE)
          snperr.pos=rep(0,nsnperr);snperr.chr=rep("",nsnperr)
          for(j in 1:nchr){
            tmp.idx=snperr.coord>=chr.start[j] & snperr.coord<=chr.end[j]
            if(sum(tmp.idx)>0){
              snperr.chr[tmp.idx]=chr.names[j]            
              snperr.pos[tmp.idx]=snperr.coord[tmp.idx]-chr.start[j]
            }
          }
          snpdet.err=data.frame(Chromosome=snperr.chr,Position=snperr.pos,
                                RefAllele=rep("error",nsnperr),AltAllele=rep(NA,nsnperr))
          out@snp.info=rbind(out@snp.info,snpdet.err)
          tmp.ord=order(out@snp.info$Chromosome,out@snp.info$Position)
          out@snp.info=out@snp.info[tmp.ord,]
          out@refallele.readcount=rbind(out@refallele.readcount,tmp[,1:npop])[tmp.ord,]
          tmp=tmp[,-1*(1:npop)]
          out@readcoverage=rbind(out@readcoverage,tmp)[tmp.ord,]
          rm(tmp) ; gc()
        }
      }
  }
  out@nsnp=nrow(out@snp.info)
  if(verbose){cat("Simulation finished:",out@nsnp,"retained in the output pooldata object\n")}
  return(out)
}


