#' Compute pairwise population population FST matrix (and possibly all pairwise SNP-specific FST)
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param method Either "Anova" (default method as described in the manuscript) or "Identity" (relies on an alternative modeling consisting in estimating unbiased Probability of Identity within and across pairs of pools)
#' @param min.cov.per.pool For Pool-Seq data (i.e., pooldata objects) only: minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons
#' @param max.cov.per.pool For Pool-Seq data (i.e., pooldata objects) only: maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons.
#' @param min.indgeno.per.pop  For allele count data (i.e., countdata objects) only: minimal number of overall counts required in each population. If at least one pop is not genotyped for at least min.indgeno.per.pop (haploid) individual, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio over all read counts for the reference allele over the read coverage)  in the pairwise comparisons.
#' @param output.snp.values If TRUE, provide SNP-specific pairwise FST for each comparisons (may lead to a huge result object if the number of pools and/or SNPs is large)
#' @param nsnp.per.bjack.block Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling) 
#' @param verbose If TRUE extra information is printed on the terminal
#' @return An object of class pairwisefst (see help(pairwisefst) for details)
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  PairwiseFST=compute.pairwiseFST(pooldata)
#' @export
compute.pairwiseFST<-function(x,method="Anova",min.cov.per.pool=-1,max.cov.per.pool=1e6,min.indgeno.per.pop=-1,min.maf=-1,output.snp.values=FALSE,nsnp.per.bjack.block=0,verbose=TRUE){
  if(!(method %in% c("Identity","Anova"))){stop("method should either be Identity or Anova (default)")}
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")
  } 
  if(is.pooldata(x)){npops=x@npools ; popnames=x@poolnames ; pooldata=TRUE}
  if(is.countdata(x)){npops=x@npops ; popnames=x@popnames ; pooldata=FALSE}  

  if(npops<3){stop("At least 3 pools (or pops) should be given to compute the matrix")}
  if(npops>50){cat("WARNING: with large number of pop samples you may also consider using the compute.fstats function (significantly faster and far more memory-efficient) to obtain the pairwise FST matrix (pairwise.fst slot of the fstats output object) that can be visualized using heatmap or other conventional clustering techniques (see the vignette).\n")}
  time1=proc.time()
    
  out=new("pairwisefst")
  out@PairwiseFSTmatrix=matrix(NA,npops,npops)
  colnames(out@PairwiseFSTmatrix)=rownames(out@PairwiseFSTmatrix)=popnames  
  n.pairs=npops*(npops-1)/2
  pop.pairs=t(combn(npops,2)) #may be usefull is output.snp.values or jacknife
  pop.pairs.names=paste0(popnames[pop.pairs[,1]],";",popnames[pop.pairs[,2]])
  out@values=as.data.frame(matrix(0,n.pairs,7))
  colnames(out@values)=c("Fst Estimate","Fst bjack mean","Fst bjack s.e.","Q2 Estimate","Q2 bjack mean","Q2 bjack s.e.","Nsnp")
  rownames(out@values)=pop.pairs.names
  out@blockjacknife=ifelse(nsnp.per.bjack.block>0,TRUE,FALSE)

  if(output.snp.values){
    tmp.size.GB=x@nsnp*n.pairs*8/1e9
    if(tmp.size.GB>4){
      cat("WARNING: because output.snp.values=TRUE, three large data frame each of >4 GB will be created to store SNP-specific Q1, Q2 and Fst for each population pairs. Be sure to have enough RAM available and that you really need SNP-specific estimates (that also can be obtained with computeFST function for a given specific pair or over all populations)\n")
    }
    out@PairwiseSnpQ1=matrix(NA,x@nsnp,n.pairs)#,dimnames=list(paste0("rs",1:x@nsnp),pop.pairs.names))
    colnames(out@PairwiseSnpQ1)=pop.pairs.names
    out@PairwiseSnpQ2=out@PairwiseSnpQ1
    out@PairwiseSnpFST=out@PairwiseSnpQ1
  }
  #crit cov snp per pop
  if(pooldata){
    snp.cov.filt=x@readcoverage>=min.cov.per.pool & x@readcoverage<=max.cov.per.pool
  }else{
    snp.cov.filt=x@total.count>=min.indgeno.per.pop
  }
  ###
  if(verbose){
    cat("Computation of the",n.pairs,"pairwise Fst\n")
    pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.pairs, clear = FALSE, width= 60)
  }
 ###
  cnt.pair=0
  pair.code=matrix(c(0,1),nrow=1) #for built-in cpp identity function
  for(p in 1:(npops-1)){
    for(q in (p+1):(npops)){
     cnt.pair=cnt.pair+1
     ##snp.sel
     if(pooldata){
       tmp.Ntot=rowSums(x@readcoverage[,c(p,q)])
       tmp.maf=0.5-abs(0.5-rowSums(x@refallele.readcount[,c(p,q)])/tmp.Ntot)
     }else{
       tmp.Ntot=rowSums(x@total.count[,c(p,q)])
       tmp.maf=0.5-abs(0.5-rowSums(x@refallele.count[,c(p,q)])/tmp.Ntot)
     }
     tmp.snp.sel=tmp.Ntot>0 & snp.cov.filt[,p] & snp.cov.filt[,q] & (tmp.maf>min.maf) 
     tmp.nsnp=sum(tmp.snp.sel)
     ##calcul snp.Q1 et snp.Q2 selon methode et type de donnees
     if(method=="Identity"){
       if(is.countdata(x)){
         tmp.Q1=.compute_snpQ1(x@refallele.count[tmp.snp.sel,c(p,q)],x@total.count[tmp.snp.sel,c(p,q)],rep(1,2),verbose=FALSE)
         tmp.Q2=.compute_snpQ2(x@refallele.count[tmp.snp.sel,c(p,q)],x@total.count[tmp.snp.sel,c(p,q)],pair.code,verbose=FALSE)
       }else{
         tmp.Q1=.compute_snpQ1(x@refallele.readcount[tmp.snp.sel,c(p,q)],x@readcoverage[tmp.snp.sel,c(p,q)],
                               x@poolsizes[c(p,q)]/(x@poolsizes[c(p,q)]-1),verbose=FALSE) 
         tmp.Q2=.compute_snpQ2(x@refallele.readcount[tmp.snp.sel,c(p,q)],x@readcoverage[tmp.snp.sel,c(p,q)],pair.code,verbose=FALSE)
       }
       num.fst=tmp.Q1 - tmp.Q2
       den.fst=1. - tmp.Q2 #same as den.fgt
        }

     if(method=="Anova"){
       if(is.countdata(x)){
         tmp.aov=.compute_snpFstAov(x@refallele.count[tmp.snp.sel,c(p,q)],x@total.count[tmp.snp.sel,c(p,q)],rep(-1,2),verbose=FALSE)
       }else{
         tmp.aov=.compute_snpFstAov(x@refallele.readcount[tmp.snp.sel,c(p,q)],x@readcoverage[tmp.snp.sel,c(p,q)],
                                    x@poolsizes[c(p,q)],verbose=FALSE)
       }
       tmp.Q1=1.-tmp.aov[,1]                                         #1-MSG
       tmp.Q2=tmp.Q1 - (tmp.aov[,2]-tmp.aov[,1])/tmp.aov[,3] #1 - MSG - (MSP - MSG) / n_c       
       num.fst=tmp.aov[,2]-tmp.aov[,1]                          #MSP-MSG :
       den.fst=tmp.aov[,2]+(tmp.aov[,3]-1)*tmp.aov[,1]          #MSP + (n_c- 1) * MSG 
       #important de ne pas passer par les Q1 et Q2 pour garder le facteur nc le calcul multilocus      
     }
     ## 
     keep=!(is.na(num.fst) | is.na(den.fst)) 
     out@values[cnt.pair,1]=out@PairwiseFSTmatrix[p,q]=out@PairwiseFSTmatrix[q,p]=sum(num.fst[keep])/sum(den.fst[keep])
     out@values[cnt.pair,4]=mean(tmp.Q2[keep])
     out@values[cnt.pair,7]=tmp.nsnp
     if(out@blockjacknife){
       fake.obj=new("countdata")
       fake.obj@nsnp=tmp.nsnp ; fake.obj@snp.info=x@snp.info[tmp.snp.sel,] 
       bjack.blocks=generate.jackknife.blocks(fake.obj,nsnp.per.bjack.block,verbose=FALSE)
       bjack.fact=(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
       
       keep=keep & !is.na(bjack.blocks$snp.block.id)
       tmp.snp.block.id=bjack.blocks$snp.block.id[keep]-1 #0 indexed for cpp
       num.fst.blk.contrib=.block_sum(num.fst[keep],tmp.snp.block.id)
       den.fst.blk.contrib=.block_sum(den.fst[keep],tmp.snp.block.id)
       sampled.fst=(num.fst.blk.contrib-sum(num.fst.blk.contrib))/(den.fst.blk.contrib-sum(den.fst.blk.contrib))
       out@values[cnt.pair,2]=mean(sampled.fst) ; out@values[cnt.pair,3]=sd(sampled.fst)*bjack.fact
       sampled.q2=.block_sum(tmp.Q2[keep],tmp.snp.block.id)
       sampled.q2=(sum(sampled.q2)-sampled.q2)/(tmp.nsnp-nsnp.per.bjack.block)
       out@values[cnt.pair,5]=mean(sampled.q2) ; out@values[cnt.pair,6]=sd(sampled.q2)*bjack.fact
       }
     if(output.snp.values){
       out@PairwiseSnpQ1[tmp.snp.sel,cnt.pair]=tmp.Q1
       out@PairwiseSnpQ2[tmp.snp.sel,cnt.pair]=tmp.Q2
       out@PairwiseSnpFST[tmp.snp.sel,cnt.pair]=num.fst/den.fst
      }
    ###     
      if(verbose){pb$tick()}#        setTxtProgressBar(pb, cnt.pair)}
    }
  }
  if(!out@blockjacknife){out@values=out@values[,c(1,4,7)]}
  time.elapsed=((proc.time()-time1)[3])
  nhours=floor(time.elapsed/3600) ; nminutes=floor((time.elapsed-nhours*3600)/60) ;  nseconds=round(time.elapsed-nhours*3600-nminutes*60)    
  cat("\nOverall Analysis Time:",nhours,"h",nminutes, "m",nseconds,"s\n")  
  
  return(out)
}
