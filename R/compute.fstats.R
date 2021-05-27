#' Estimate the F-statistics (F2, F3, F3star, F4, Dstat)
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param nsnp.per.bjack.block Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling) 
#' @param computeDstat If TRUE compute Dstatistics (i.e. scaled F4). This may add some non negligible computation time if the number of population is large (n>15)
#' @param return.F4.blockjackknife.samples If TRUE (and nsnp.per.bjack.block>0) return F4 estimates for each block-jackknife sample (useful to compute F4 ratios standard errors)
#' @param verbose If TRUE extra information is printed on the terminal
# #' @param nocpp If  TRUE no cpp
#' @details The function estimates for the n populations (or pools) represented in the input object x:
#' \enumerate{
#' \item The F2 statistics for all the \eqn{n(n-1)/2} pairs of populations (or pools) and their scaled version (equivalent to Fst as compute with \code{\link{compute.pairwiseFST}} with method="Identity")
#' \item If n>2, The F3 statistics for all the \eqn{npools(npools-1)(npools-2)/2} possible triplets of populations (or pools) and their scaled version (named F3star after Patterson et al., 2012)
#' \item If n>3, The F4 statistics and the D-statistics (a scaled version of the F4) for all the \eqn{npools(npools-1)(npools-2)*(npools-3)/8} possible quadruplets of populations
#' \item The estimated within population heterozygosities (=1-Q1)
#' }
#' @return An object of class fstats (see help(fstats) for details)
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fstats=compute.fstats(pooldata)
#' @export
compute.fstats<-function(x,nsnp.per.bjack.block=0,computeDstat=FALSE,return.F4.blockjackknife.samples=FALSE,verbose=TRUE){
  time1=proc.time()
  nocpp=FALSE #to compare R version and cpp of some functions
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")
    } 
  if(is.pooldata(x)){npops=x@npools ; popnames=x@poolnames}
  if(is.countdata(x)){npops=x@npops ; popnames=x@popnames}  
  if(npops<4){computeF4=FALSE;computeDstat=FALSE}else{computeF4=TRUE}
  if(npops<3){computeF3=FALSE}else{computeF3=TRUE}
  
  out=new("fstats")
  if(nsnp.per.bjack.block>0){out@blockjacknife=TRUE}else{out@blockjacknife=FALSE}
  
  #Computing Q1 for each pop  
  if(verbose){cat("Estimating Q1\n")}
  if(is.pooldata(x)){
  # R1=x@refallele.readcount*(x@refallele.readcount-1)/(x@readcoverage*(x@readcoverage-1))
  # R1 = (1/(matrix(1,x@nsnp,x@npools) %*% diag(x@poolsizes-1)))*(R1 %*% diag(x@poolsizes) - x@refallele.readcount/x@readcoverage) #estimateur sans bias de Y(Y-1)/(N*(N-1))
   x.count.alt=x@readcoverage-x@refallele.readcount 
   # Q1=(x@refallele.readcount*(x@refallele.readcount-1) + x.count.alt*(x.count.alt-1))/(x@readcoverage*(x@readcoverage-1))
   # Q1 =t((t(Q1)*x@poolsizes-1)/(x@poolsizes-1)) 
   Q1=2*x@refallele.readcount*x.count.alt/(x@readcoverage*(x@readcoverage-1))
   Q1 =1 - t(t(Q1)*x@poolsizes/(x@poolsizes-1)) 
  }else{
  # R1=x@refallele.count*(x@refallele.count-1) /(x@total.count*(x@total.count-1)) 
    x.count.alt=x@total.count-x@refallele.count  
    #Q1=(x@refallele.count*(x@refallele.count-1) + x.count.alt*(x.count.alt-1))/(x@total.count*(x@total.count-1))
    Q1=1-2*x@refallele.count*x.count.alt/(x@total.count*(x@total.count-1))
  }
  
  #Computing Q2 for each snp/pop pair
  n.f2=npops*(npops-1)/2
  if(verbose){
    cat("Estimating Q2\n")
    pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.f2, clear = FALSE, width= 60)
    }
  Q2=matrix(0,x@nsnp,n.f2)
  p1p2.idx=matrix(0,n.f2,2)
  f2.idx=matrix(NA,npops,npops) #useful to store indexes of pop and idx of F2
  f2.names=rep("",n.f2)
  cnt=0
  for(i in 1:(npops-1)){
    for(j in (i+1):npops){
        cnt=cnt+1
        p1p2.idx[cnt,]=c(i,j) ; f2.idx[i,j]=f2.idx[j,i]=cnt
        f2.names[cnt]=paste(popnames[i],popnames[j],sep=",")
        Q2[,cnt]=x.count.alt[,i]*x.count.alt[,j] #x.count.alt a ete defini plus haut (differe selon objet)
        if(is.pooldata(x)){
    #     R2[,cnt]=x@refallele.readcount[,i]*x@refallele.readcount[,j]/(x@readcoverage[,i]*x@readcoverage[,j])
          Q2[,cnt]=(x@refallele.readcount[,i]*x@refallele.readcount[,j] + Q2[,cnt])/(x@readcoverage[,i]*x@readcoverage[,j]) 
        }else{
    #     R2[,cnt]=x@refallele.count[,i]*x@refallele.count[,j]/(x@total.count[,i]*x@total.count[,j]) 
          Q2[,cnt]=(x@refallele.count[,i]*x@refallele.count[,j] + Q2[,cnt])/(x@total.count[,i]*x@total.count[,j]) 
        }
        if(verbose){pb$tick()}  
      }
  }
  rm(x.count.alt) ; gc()
  if(verbose){pb$terminate()}
  
 ###Computing heterozygosity
  if(verbose){cat("Estimating within-population heterozygosities\n")}
  het.est=1-colMeans(Q1,na.rm=T)
 ####Computing F2
  if(verbose){cat("Estimating F2\n")}
  hat.Q1=colMeans(Q1,na.rm=T)
  hat.Q2=colMeans(Q2,na.rm=T)
  f2.est=(hat.Q1[p1p2.idx[,1]]+hat.Q1[p1p2.idx[,2]])/2 - hat.Q2
  fst.est=f2.est/(1-hat.Q2)
  
  ###Compute F3 from F2: F3(A;B,C)=(F2(A,B)+F2(A,C)-F2(B,C))/2
  if(computeF3){
    if(verbose){cat("Estimating F3\n")}
    n.f3=npops*(npops-1)*(npops-2)/2 #npop*choose(npop-1,2)
    f3.est=rep(0,n.f3)
    f3.names=rep("",n.f3)
    f3.f2idx=matrix(0,n.f3,3) #useful if block-jackknife
    f3.q1idx=rep(0,n.f3) #useful for F3star
    cnt=0
    for(i in 1:npops){
      for(j in 1:(npops-1)){
        for(k in (j+1):npops){
          if(j!=i & k!=i){
            cnt=cnt+1
            f3.names[cnt]=paste0(popnames[i],";",popnames[j],",",popnames[k])
            f3.f2idx[cnt,]=c(f2.idx[i,j],f2.idx[i,k],f2.idx[j,k])
            f3.q1idx[cnt]=i
          }
        }
      }
    }
   f3.est=(f2.est[f3.f2idx[,1]] + f2.est[f3.f2idx[,2]] -f2.est[f3.f2idx[,3]] )/2
   f3star.est=f3.est/(1-hat.Q1[f3.q1idx])
  }
  
  ###compute F4 from F2: F4(A,B;C,D)=(F2(A,D)+F2(B,C)-F2(A,C)-F2(B,D))/2
  if(computeF4){
    if(verbose){cat("Estimating F4\n")}
    n.f4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 #choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
    f4.est=rep(0,n.f4)
    f4.names=rep("",n.f4)
    f4.f2idx=matrix(0,n.f4,6) #useful if block-jackknife
    cnt=0
    for(i in 1:(npops-1)){
      for(j in (i+1):npops){
        tmp.popindex=(1:npops)[-c(i,j)] 
        tmp.popindex=tmp.popindex[tmp.popindex>i] 
        tmp.popindex.lg=length(tmp.popindex)
        if(tmp.popindex.lg>1){
          for(k in 1:(tmp.popindex.lg-1)){
            for(l in (k+1):tmp.popindex.lg){  
              p=tmp.popindex[k] ; q=tmp.popindex[l] 
              cnt=cnt+1 
              f4.names[cnt]=paste0(popnames[i],",",popnames[j],";",popnames[p],",",popnames[q])
              f4.f2idx[cnt,]=c(f2.idx[i,q],f2.idx[j,p],f2.idx[i,p],f2.idx[j,q],f2.idx[i,j],f2.idx[p,q])
            }}}}}  
    f4.est=(f2.est[f4.f2idx[,1]] + f2.est[f4.f2idx[,2]] - f2.est[f4.f2idx[,3]] - f2.est[f4.f2idx[,4]])/2
    if(computeDstat){
      if(verbose){cat(" Computing Dstat\n")}
      if(nocpp){
#    d.denom=colMeans((1-Q2[,f4.f2idx[,5]])*(1-Q2[,f4.f2idx[,6]]),na.rm=T)
     if(verbose){pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.f4, clear = FALSE, width= 60)}
    d.denom=c() #oblige de faire une boucle car consomme trop de memoire sinon
    for(i in 1:n.f4){
      d.denom=c(d.denom,mean((1-Q2[,f4.f2idx[i,5]])*(1-Q2[,f4.f2idx[i,6]]),na.rm=T))
      if(verbose){pb$tick()}
    }
    if(verbose){pb$terminate()}
      }else{
      d.denom=.compute_Ddenom(Q2,f4.f2idx[,5:6],verbose)
      }
      d.est=f4.est/d.denom
    }
  }
 
  #######Jackknife
  ##By definition: se(hat(f))={Sum((fi-mean_f)^2)*(n-1)/n}^0.5={var_f*(n-1)*(n-1)/n}^0.5=sd_f*(n-1)/(n^0.5)
  if(out@blockjacknife){
    if(verbose){cat("Starting Block-Jackknife sampling\n")}
    bjack.blocks=generate.jackknife.blocks(x,nsnp.per.bjack.block,verbose=verbose)
    tmp.sel=!is.na(bjack.blocks$snp.block.id)
    Q1=Q1[tmp.sel,] ; gc()
    Q2=Q2[tmp.sel,] ; gc()
    block.id=bjack.blocks$snp.block.id[tmp.sel]

    if(verbose){cat("   computing Q1 averages per blocks\n")} 

    if(nocpp){    
    tmp.Q1.per.block=apply(Q1,2,function(y){return(by(y,block.id,sum,na.rm=T))})
    tmp.Q1.per.block=colSums(tmp.Q1.per.block,na.rm=T)-t(tmp.Q1.per.block)
    tmp.nsnp.per.sample=apply(1-is.na(Q1),2,function(y){return(by(y,block.id,sum))})
    tmp.nsnp.per.sample=colSums(tmp.nsnp.per.sample,na.rm=T)-t(tmp.nsnp.per.sample)
    sampled.Q1=tmp.Q1.per.block/tmp.nsnp.per.sample
    }else{
    sampled.Q1=.compute_Q_bjmeans(Q1,block.id,verbose)
    }
    
    if(verbose){cat("   computing Q2 averages per blocks\n")}    

    if(nocpp){      
    tmp.Q2.per.block=apply(Q2,2,function(y){return(by(y,block.id,sum,na.rm=T))})
    tmp.Q2.per.block=colSums(tmp.Q2.per.block,na.rm=T)-t(tmp.Q2.per.block)
    tmp.nsnp.per.sample=apply(1-is.na(Q2),2,function(y){return(by(y,block.id,sum))})
    tmp.nsnp.per.sample=colSums(tmp.nsnp.per.sample,na.rm=T)-t(tmp.nsnp.per.sample)
    sampled.Q2=tmp.Q2.per.block/tmp.nsnp.per.sample
    }else{
    sampled.Q2=.compute_Q_bjmeans(Q2,block.id,verbose)
    }
    
    if(verbose){cat("   computing F2 averages per blocks\n")}    

    if(nocpp){      
    F2=(Q1[,p1p2.idx[,1]]+Q1[,p1p2.idx[,2]])/2 - Q2
    tmp.nsnp.per.sample=apply(1-is.na(F2),2,function(y){return(by(y,block.id,sum))})
    tmp.nsnp.per.sample=colSums(tmp.nsnp.per.sample,na.rm=T)-t(tmp.nsnp.per.sample)
    rm(Q1) ; gc()
    tmp.F2.per.block=apply(F2,2,function(y){return(by(y,block.id,sum,na.rm=T))})
    rm(F2) ; gc()
    tmp.F2.per.block=colSums(tmp.F2.per.block,na.rm=T)-t(tmp.F2.per.block)
    sampled.f2=tmp.F2.per.block/tmp.nsnp.per.sample
    }else{   
    sampled.f2=.compute_F2_bjmeans(Q1,Q2,p1p2.idx,block.id,verbose)
    rm(Q1) ; gc()
    if(!computeDstat){rm(Q2) ; gc()}
    }
    
   # if(computeDstat){        
   #  if(verbose){cat("   computing Dstat denominator averages per blocks\n")}
   #   if(nocpp){
   #       if(verbose){pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.f4, clear = FALSE, width= 60)}
   #  #par securite on fait une boucle car D.denom peut devenir tres gros (et explose la memoire) si beaucoup de pops
   #  sampled.D.denom=matrix(0,n.f4,ncol(sampled.f2))
   #  for(i in 1:n.f4){
   #    tmp.D.denom=(1-Q2[,f4.f2idx[i,5]])*(1-Q2[,f4.f2idx[i,6]])
   #    tmp.nsnp.per.sample=by(1-is.na(tmp.D.denom),block.id,sum)
   #    tmp.nsnp.per.sample=sum(tmp.nsnp.per.sample)-tmp.nsnp.per.sample
   #    tmp.D.denom.per.block=by(tmp.D.denom,block.id,sum,na.rm=T)
   #    tmp.D.denom.per.block=sum(tmp.D.denom.per.block)-tmp.D.denom.per.block
   #    sampled.D.denom[i,]=tmp.D.denom.per.block/tmp.nsnp.per.sample
   #    if(verbose){pb$tick()}
   #  }
   #  if(verbose){pb$terminate()}
   #   }else{
   #   sampled.D.denom=.compute_Ddenom_bjmeans(Q2,f4.f2idx[,5:6],block.id,verbose)
   #   }
   # }
   #  rm(Q2) ; gc()

    if(verbose){cat("Starting computation of estimators s.e.\n")}       
        
    ##compute heterozygosity from Q1
    sampled.het=1 - sampled.Q1
    out@heterozygosities=as.data.frame(cbind(het.est,rowMeans(sampled.het),apply(sampled.het,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)))
    rownames(out@heterozygosities)=popnames
    colnames(out@heterozygosities)=c("Estimate","bjack mean","bjack s.e.")
    rm(sampled.het) ; gc() 

    if(verbose){cat("   within-pop heterozygosity s.e. estimation done\n")}       
        
    ##F2 from Q1 and Q2
    out@f2.values=as.data.frame(cbind(f2.est,rowMeans(sampled.f2),apply(sampled.f2,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)))
    rownames(out@f2.values)=f2.names
    colnames(out@f2.values)=c("Estimate","bjack mean","bjack s.e.")
   #fst
    sampled.fst=sampled.f2/(1-sampled.Q2)
    out@fst.values=as.data.frame(cbind(fst.est,rowMeans(sampled.fst),apply(sampled.fst,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)))
    rownames(out@fst.values)=f2.names
    colnames(out@fst.values)=c("Estimate","bjack mean","bjack s.e.")
    rm(sampled.fst) ; gc()

    if(verbose){cat("   F2 s.e. estimation done\n")}     
        
  if(computeF3){
     tmp.f3.est=tmp.f3star.est=matrix(0,n.f3,2)
     sampled.f3=(sampled.f2[f3.f2idx[,1],] + sampled.f2[f3.f2idx[,2],] -sampled.f2[f3.f2idx[,3],] )/2
     tmp.f3.est[,1]=rowMeans(sampled.f3)
     tmp.f3.est[,2]=apply(sampled.f3,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
     #computation of Q.matrix
     tmp=rbind(sampled.f2,sampled.f3)
     try(out@Q.matrix<-cov(t(tmp))*((bjack.blocks$nblocks-1)**2)/bjack.blocks$nblocks,silent=TRUE)
     if(nrow(out@Q.matrix)==0){
       cat("Warning: not enough memory to compute and store the error covariance matrix (Q.matrix).\nIt will not be possible to fit graph with the output fstats object\nYou may try to run compute.fstats on a smaller object, i.e., with less populations\n(use pooldata.subset or countdata.subset to create it) involving only the populations of interest.\n")
     }else{
      rownames(out@Q.matrix)=colnames(out@Q.matrix)=c(f2.names,f3.names)
      rm(tmp) ; gc()
     }
     #F3 star
     sampled.f3=sampled.f3/(1-sampled.Q1[f3.q1idx,])
     tmp.f3star.est[,1]=rowMeans(sampled.f3)
     tmp.f3star.est[,2]=apply(sampled.f3,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)     
     rm(sampled.f3) ; gc()
     #filling out object
     out@f3.values=as.data.frame(cbind(f3.est,tmp.f3.est,tmp.f3.est[,1]/tmp.f3.est[,2]))
     colnames(out@f3.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")
     rownames(out@f3.values)=f3.names

     out@f3star.values=as.data.frame(cbind(f3star.est,tmp.f3star.est,tmp.f3star.est[,1]/tmp.f3star.est[,2]))
     colnames(out@f3star.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")
     rownames(out@f3star.values)=f3.names

     rm(tmp.f3.est,tmp.f3star.est) ; gc()
     if(verbose){cat("   F3 and F3* s.e. estimation done\n")} 
  }

  if(computeF4){
    if(verbose){
      if(computeDstat){
        cat("   estimating F4 and Dstat s.e. (may be long since require denominator averages per blocks)\n")
      }else{cat("   estimating F4 s.e.\n")}
      }
    tmp.f4.est=matrix(0,n.f4,2)
    if(computeDstat){tmp.D.est=matrix(0,n.f4,2)}
      #on segmente par blocks de npop*(npop-1)/2 pops : de npops=5 a 500: les blocs ont toujours plus de 10 (pas de pb avecmatrix) 
      if(npops>10){
       id.begin=seq(1,n.f4,npops*(npops-1)/2)
       id.end=c(id.begin[-1]+1,n.f4)
      }else{
        id.begin=1 ; id.end=n.f4
      }
      if(verbose){
        pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = length(id.begin), clear = FALSE, width= 60)
      }
      for(i in 1:length(id.begin)){
        tmp.idx=id.begin[i]:id.end[i]        
        sampled.f4<-(sampled.f2[f4.f2idx[tmp.idx,1],] + sampled.f2[f4.f2idx[tmp.idx,2],] - sampled.f2[f4.f2idx[tmp.idx,3],] - sampled.f2[f4.f2idx[tmp.idx,4],])/2
        tmp.f4.est[tmp.idx,1]=rowMeans(sampled.f4)
        tmp.f4.est[tmp.idx,2]=apply(sampled.f4,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
        if(computeDstat){
          sampled.f4=sampled.f4/.compute_Ddenom_bjmeans(Q2,f4.f2idx[tmp.idx,5:6],block.id,FALSE)
          tmp.D.est[tmp.idx,1]=rowMeans(sampled.f4)
          tmp.D.est[tmp.idx,2]=apply(sampled.f4,1,sd)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
        }
        if(verbose){pb$tick()}
      }
      if(computeDstat){rm(Q2);gc()}
      rm(sampled.f4) ; gc()
      if(verbose){pb$terminate()}
      ##out objects
      out@f4.values=as.data.frame(cbind(f4.est,tmp.f4.est,tmp.f4.est[,1]/tmp.f4.est[,2]))
      colnames(out@f4.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")
      rownames(out@f4.values)=f4.names
      rm(tmp.f4.est) ; gc()
      if(computeDstat){
        out@Dstat.values=as.data.frame(cbind(d.est,tmp.D.est,tmp.D.est[,1]/tmp.D.est[,2]))
        colnames(out@Dstat.values)=c("Estimate","bjack mean","bjack s.e.","Z-score")
        rownames(out@Dstat.values)=f4.names
        rm(tmp.D.est) ; gc()
        if(verbose){cat("   F4 and D s.e. estimation done\n")} 
      }else{
        if(verbose){cat("   F4 s.e. estimation done\n")} 
      }
    ###samples pour F4 ratio (test si faisbale)  
    ##stockage: a la fin car sinon bouffe de la memoire necessaire pendant le calcul (recalcul pas trop couteux en fait)
      if(return.F4.blockjackknife.samples){
        try(out@F4.bjack.samples<- (sampled.f2[f4.f2idx[,1],] + sampled.f2[f4.f2idx[,2],] - sampled.f2[f4.f2idx[,3],] - sampled.f2[f4.f2idx[,4],])/2,silent=TRUE)
        if(nrow(out@F4.bjack.samples)==0){#test dimensions        
          cat("Warning: not enough memory to store the F4 block-jackknife samples: options disregarded.\nIf really needed to compute F4 ratio, you may try to run compute.fstats on a smaller object (use pooldata.subset or countdata.subset to create it) involving only the populations of interest.\n")
          return.F4.blockjackknife.samples=FALSE
        }else{
          rownames(out@F4.bjack.samples)=f4.names 
        }
      }  
  }
    
  }else{
    out@f2.values=as.data.frame(matrix(f2.est,ncol=1,dimnames = list(f2.names,"Estimate")))
    out@fst.values=as.data.frame(matrix(fst.est,ncol=1,dimnames = list(f2.names,"Estimate")))    
    out@heterozygosities=as.data.frame(matrix(het.est,ncol=1,dimnames = list(popnames,"Estimate")))
    if(computeF3){
      out@f3.values=as.data.frame(matrix(f3.est,ncol=1,dimnames = list(f3.names,"Estimate")))
      out@f3star.values=as.data.frame(matrix(f3star.est,ncol=1,dimnames = list(f3.names,"Estimate")))      
      }
    if(computeF4){
      out@f4.values=as.data.frame(matrix(f4.est,ncol=1,dimnames = list(f4.names,"Estimate")))
      if(computeDstat){out@Dstat.values=as.data.frame(matrix(d.est,ncol=1,dimnames = list(f4.names,"Estimate")))}      
      }
  }
  
  ##names of the comparison
  out@comparisons[["pops"]]=popnames
  out@comparisons[["F2"]]=out@comparisons[["Fst"]]=matrix(unlist(strsplit(rownames(out@f2.values),split="[,;]")),nrow=nrow(out@f2.values),byrow=T,dimnames=list(rownames(out@f2.values),c("P1","P2")))
  if(computeF3){
    out@comparisons[["F3"]]=out@comparisons[["F3star"]]=matrix(unlist(strsplit(rownames(out@f3.values),split="[,;]")),nrow=nrow(out@f3.values),byrow=T,dimnames=list(rownames(out@f3.values),c("Px","P1","P2")))
    }
  if(computeF4){
    out@comparisons[["F4"]]=matrix(unlist(strsplit(rownames(out@f4.values),split="[,;]")),nrow=nrow(out@f4.values),byrow=T,dimnames=list(rownames(out@f4.values),c("P1","P2","P3","P4")))}  
  if(computeDstat){out@comparisons[["Dstat"]]=out@comparisons[["F4"]]}

  if(verbose){
  time.elapsed=((proc.time()-time1)[3])
  nhours=floor(time.elapsed/3600) ; nminutes=floor((time.elapsed-nhours*3600)/60) ;  nseconds=round(time.elapsed-nhours*3600-nminutes*60)  
  cat("\nOverall Analysis Time:",nhours,"h",nminutes, "m",nseconds,"s\n")  
  }
  
  return(out)
  
}
