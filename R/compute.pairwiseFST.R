#' Compute pairwise population population FST matrix (and possibly all pairwise SNP-specific FST)
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param method Either "Anova" (default method as described in the manuscript) or "Identity" (relies on an alternative modeling consisting in estimating unbiased Probability of Identity within and across pairs of pools)
#' @param min.cov.per.pool For Pool-Seq data (i.e., pooldata objects) only: minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons
#' @param max.cov.per.pool For Pool-Seq data (i.e., pooldata objects) only: maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded in the corresponding pairwise comparisons.
#' @param min.indgeno.per.pop  For allele count data (i.e., countdata objects) only: minimal number of overall counts required in each population. If at least one pop is not genotyped for at least min.indgeno.per.pop (haploid) individual, the position is discarded
#' @param min.maf Minimal allowed Minor Allele Frequency (computed from the ratio overal read counts for the reference allele over the read coverage)  in the pairwise comparisons.
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
  if(npops>20){warning("with >20 pop samples it is recommended to utilize the compute.fstats function (significantly faster and far more memory-efficient) to obtain the pairwise FST matrix (pairwise.fst slot of the fstats output object) that can be visualized using heatmap or other conventional clustering techniques (see the vignette).")}
  time1=proc.time()
    
  out=new("pairwisefst")
  if(nsnp.per.bjack.block>0){out@blockjacknife=TRUE}else{out@blockjacknife=FALSE}
  out@PairwiseFSTmatrix=matrix(NA,npops,npops)
  colnames(out@PairwiseFSTmatrix)=rownames(out@PairwiseFSTmatrix)=popnames  
  n.pairs=npops*(npops-1)/2
  pop.pairs=t(combn(npops,2)) #may be usefull is output.snp.values or jacknife
  pop.pairs.names=paste0(popnames[pop.pairs[,1]],";",popnames[pop.pairs[,2]])
  out@values=as.data.frame(matrix(0,n.pairs,7))
  colnames(out@values)=c("Fst Estimate","Fst bjack mean","Fst bjack s.e.","Q2 Estimate","Q2 bjack mean","Q2 bjack s.e.","Nsnp")
  rownames(out@values)=pop.pairs.names
  if(output.snp.values | out@blockjacknife){
    out@PairwiseSnpQ1=matrix(NA,x@nsnp,n.pairs)#,dimnames=list(paste0("rs",1:x@nsnp),pop.pairs.names))
    colnames(out@PairwiseSnpQ1)=pop.pairs.names
    out@PairwiseSnpQ2=out@PairwiseSnpQ1
    if(output.snp.values){out@PairwiseSnpFST=out@PairwiseSnpQ1}
    if(method=="Anova"){snp.Nc=matrix(NA,x@nsnp,n.pairs)}
    #Si anova: FST=sum(Nc*(snpQ1-snpQ2))/sum(Nc*(1-snpQ2))=sum(MSP-MSI)/sum(MSP-(nc-1)MSI) #En effet nc varie entre les SNPs
  }
  ###
  #On recalcule les Q1 et Q2 pour chaque paire
  cnt.pair=0 
  if(verbose){
    cat("Computation of the",n.pairs,"pairwise Fst\n")
 #   pb <- txtProgressBar(min = 0, max = n.pairs, style = 3)
    pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.pairs, clear = FALSE, width= 60)
  }
  
  #crit cov snp per pop
  if(pooldata){
    snp.cov.filt=x@readcoverage>=min.cov.per.pool & x@readcoverage<=max.cov.per.pool
  }else{
    snp.cov.filt=x@total.count>=min.indgeno.per.pop
  }
  #optimisation: stockage produit des comptages (pas les meme selon methode)
  if(method=="Identity"){
    if(pooldata){#pour Q1
      tmp.cntalt=x@readcoverage-x@refallele.readcount
      prodcount=(x@refallele.readcount*(x@refallele.readcount-1) + tmp.cntalt*(tmp.cntalt-1))/(x@readcoverage*(x@readcoverage-1))
    }else{#pour le Q1
      tmp.cntalt=x@total.count-x@refallele.count
      prodcount=(x@refallele.count*(x@refallele.count-1) + tmp.cntalt*(tmp.cntalt-1))/(x@total.count*(x@total.count-1))
    }
    rm(tmp.cntalt)
  }else{
    if(pooldata){#pour le SSI
      prodcount=x@refallele.readcount*(x@readcoverage-x@refallele.readcount)/x@readcoverage
    }else{#pour le MSG
      prodcount=x@refallele.count*(x@total.count-x@refallele.count)/x@total.count
    }
  }
      
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
      ##calcul snp.Q1 et snp.Q2 selon methode et type de donnees
      if(method=="Identity"){
        if(!pooldata){
          snp.Q1=rowMeans(prodcount[,c(p,q)],na.rm=T) 
          snp.Q2=( x@refallele.count[,p]*x@refallele.count[,q] + (x@total.count[,p]-x@refallele.count[,p])*(x@total.count[,q]-x@refallele.count[,q]))/(x@total.count[,p]*x@total.count[,q])
        }else{#see eq. A39 (Hivert et al., 2018)
          #  Q1 = (1/(matrix(1,x@nsnp,x@npools) %*% diag(x@poolsizes-1)))*(Q1 %*% diag(x@poolsizes) - 1) 
          snp.Q1=t((t(prodcount[,c(p,q)])*x@poolsizes[c(p,q)]-1)/(x@poolsizes[c(p,q)]-1)) #equivalent optimise (utilise operation matrix*vectuer qui se fait en colonnes)
          lambdaj=x@poolsizes[c(p,q)]*(x@poolsizes[c(p,q)]-1)
          lambdaj=lambdaj/sum(lambdaj)
          snp.Q1=rowSums(t(t(snp.Q1)*lambdaj))
          snp.Q2=( x@refallele.readcount[,p]*x@refallele.readcount[,q] + (x@readcoverage[,p]-x@refallele.readcount[,p])*(x@readcoverage[,q]-x@refallele.readcount[,q]))/(x@readcoverage[,p]*x@readcoverage[,q])
        }
        hat.Q2=mean(snp.Q2[tmp.snp.sel],na.rm=TRUE)
        hat.fst=sum(snp.Q1[tmp.snp.sel]-snp.Q2[tmp.snp.sel],na.rm=TRUE)/sum(1-snp.Q2[tmp.snp.sel],na.rm=TRUE)
      }
      if (method=="Anova"){
        if(!pooldata){#eq 5.2 in Weir 1996
          SumNi=rowSums(x@total.count[,c(p,q)])
          Nic=x@total.count[,c(p,q)]-(x@total.count[,c(p,q)]**2)/SumNi
          Nc=rowSums(Nic)#/(x@npops-1) : ici 2-1=1
          MSG=2*(rowSums(prodcount[,c(p,q)])) /(SumNi-2)
          PA=rowSums(x@refallele.count[,c(p,q)])/SumNi
          MSP=2*(rowSums(x@total.count[,c(p,q)]*((x@refallele.count[,c(p,q)]/x@total.count[,c(p,q)]-PA)**2)))#/(x@npops-1) : ici 2-1=1
          snp.Q1 = 1 - MSG ; snp.Q2 = 1 - MSG - (MSP - MSG) / Nc
          hat.fst=sum(Nc[tmp.snp.sel]*(snp.Q1[tmp.snp.sel]-snp.Q2[tmp.snp.sel]),na.rm=TRUE)/sum(Nc[tmp.snp.sel]*(1-snp.Q2[tmp.snp.sel]),na.rm=TRUE)
          if(out@blockjacknife){snp.Nc[tmp.snp.sel,cnt.pair] = Nc[tmp.snp.sel] }#important de laisser a NA les non sel (ponderation automatique dans block jackknife)
        }else{#Hivert et al. 2018
          mtrx.n_i <- matrix(x@poolsizes[c(p,q)],nrow = x@nsnp,ncol = 2,byrow = TRUE)
          C_1 <- rowSums(x@readcoverage[,c(p,q)])
          C_2 <- rowSums(x@readcoverage[,c(p,q)]^2)
          D_2 <- rowSums(x@readcoverage[,c(p,q)] / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i,na.rm = TRUE)
          D_2.star <- rowSums(x@readcoverage[,c(p,q)] * (x@readcoverage[,c(p,q)] / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i),na.rm = TRUE) / C_1
          n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
          SSI <- 2*rowSums(prodcount[,c(p,q)])
          SA<-rowSums(x@refallele.readcount[,c(p,q)])/C_1
          SSP<-2*rowSums(x@readcoverage[,c(p,q)]*((x@refallele.readcount[,c(p,q)]/x@readcoverage[,c(p,q)]-SA)**2))
          MSI <- SSI / (C_1 - D_2)
          MSP <- SSP / (D_2 - D_2.star)
          snp.Q1 = 1 - MSI ; snp.Q2 = 1 - MSI - (MSP - MSI) / n_c
          hat.fst=sum(n_c[tmp.snp.sel]*(snp.Q1[tmp.snp.sel]-snp.Q2[tmp.snp.sel]),na.rm=TRUE)/sum(n_c[tmp.snp.sel]*(1-snp.Q2[tmp.snp.sel]),na.rm=TRUE)
          if(out@blockjacknife){snp.Nc[tmp.snp.sel,cnt.pair] = n_c[tmp.snp.sel] }#important de laisser a NA les non sel (ponderation automatique dans block jackknife)
        }
        hat.Q2=mean(snp.Q2[tmp.snp.sel],na.rm=TRUE)
      }
      out@PairwiseFSTmatrix[p,q]=out@PairwiseFSTmatrix[q,p]=hat.fst
      out@values[cnt.pair,1]=hat.fst
      out@values[cnt.pair,4]=hat.Q2     
      out@values[cnt.pair,7]=length(tmp.snp.sel)
      if(output.snp.values | out@blockjacknife){
        out@PairwiseSnpQ1[tmp.snp.sel,cnt.pair]=snp.Q1[tmp.snp.sel]
        out@PairwiseSnpQ2[tmp.snp.sel,cnt.pair]=snp.Q2[tmp.snp.sel]
        if(output.snp.values){out@PairwiseSnpFST[tmp.snp.sel,cnt.pair]=(snp.Q1[tmp.snp.sel]-snp.Q2[tmp.snp.sel])/(1-snp.Q2[tmp.snp.sel])}
      }      
      if(verbose){pb$tick()}#        setTxtProgressBar(pb, cnt.pair)}
    }
  }
  if(verbose){pb$terminate()}#close(pb)}
  
  ##start blockjackknife
  if(out@blockjacknife){
    if(verbose){cat("Starting Block-Jackknife sampling\n")}
    bjack.blocks=generate.jackknife.blocks(x,nsnp.per.bjack.block,verbose=verbose)
    tmp.idx.sel=!is.na(bjack.blocks$snp.block.id) 
    tmp.snp.block.id=bjack.blocks$snp.block.id[tmp.idx.sel]
    if(method=="Identity"){
      tmp.sampled.q1=apply(out@PairwiseSnpQ1[tmp.idx.sel,],2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
    }else{
      tmp.sampled.q1=apply(out@PairwiseSnpQ1[tmp.idx.sel,]*snp.Nc[tmp.idx.sel,],2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})      
    }
    if(!output.snp.values){out@PairwiseSnpQ1=matrix(NA,0,0); gc()}
    tmp.sampled.q1=colSums(tmp.sampled.q1,na.rm=T)-t(tmp.sampled.q1)  
    tmp.sampled.q2=apply(out@PairwiseSnpQ2[tmp.idx.sel,],2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
    tmp.sampled.q2=colSums(tmp.sampled.q2,na.rm=T)-t(tmp.sampled.q2)  
    if(method=="Anova"){
      tmp.sampled.ncq2=apply(out@PairwiseSnpQ2[tmp.idx.sel,]*snp.Nc[tmp.idx.sel,],2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
      tmp.sampled.ncq2=colSums(tmp.sampled.ncq2,na.rm=T)-t(tmp.sampled.ncq2)
      tmp.sampled.sumnc=apply(snp.Nc[tmp.idx.sel,],2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
      tmp.sampled.sumnc=colSums(tmp.sampled.sumnc,na.rm=T)-t(tmp.sampled.sumnc) 
      sampled.fst=(tmp.sampled.q1-tmp.sampled.ncq2) / (tmp.sampled.sumnc-tmp.sampled.ncq2)  
      rm(tmp.sampled.ncq2,tmp.sampled.sumnc)
    }else{
      tmp.snp.cnt=apply(1-(is.na(out@PairwiseSnpQ2[tmp.idx.sel,])|is.na(out@PairwiseSnpQ2[tmp.idx.sel,])),2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
      tmp.snp.cnt=colSums(tmp.snp.cnt,na.rm=T)-t(tmp.snp.cnt) 
      sampled.fst=(tmp.sampled.q1-tmp.sampled.q2) / (tmp.snp.cnt-tmp.sampled.q2)  
    }
    rm(tmp.sampled.q1)
    out@values[,2]=rowMeans(sampled.fst)
    out@values[,3]=apply(sampled.fst,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
    #correction des q2: a ce stade sum: on doit faire moyenne par blocs
    tmp.snp.cnt=apply(1-is.na(out@PairwiseSnpQ2[tmp.idx.sel,]),2,function(y){return(by(y,tmp.snp.block.id,sum,na.rm=T))})
    tmp.snp.cnt=colSums(tmp.snp.cnt,na.rm=T)-t(tmp.snp.cnt) 
    if(!output.snp.values){out@PairwiseSnpQ2=matrix(NA,0,0) ; gc()}
    tmp.sampled.q2=tmp.sampled.q2/tmp.snp.cnt
    out@values[,5]=rowMeans(tmp.sampled.q2)
    out@values[,6]=apply(tmp.sampled.q2,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks) 
  }else{
    out@values=out@values[,c(1,4,7)]
  }
  ###

  time.elapsed=((proc.time()-time1)[3])
  nhours=floor(time.elapsed/3600) ; nminutes=floor((time.elapsed-nhours*3600)/60) ;  nseconds=round(time.elapsed-nhours*3600-nminutes*60)    
  cat("\nOverall Analysis Time:",nhours,"h",nminutes, "m",nseconds,"s\n")  
  
  return(out)
}