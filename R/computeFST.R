#' Compute FST from Pool-Seq data or Count data
#' @param x A pooldata object containing Pool-Seq information or countdata object containing allele counts information
#' @param method Either "Anova" (default method as described in Hivert et al (2018, eq. 9) for pool-seq data and Weir (1996, eq. 5.2) for count data) or "Identity" (relying on unbiased estimators of Probability of Identity within and across pairs of pools/populations)
#' @param nsnp.per.bjack.block Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling) 
#' @param sliding.window.size Number of consecutive SNPs within a window for multi-locus computation of Fst over sliding window with half-window size step (default=0, i.e., no sliding-window scan) 
#' @param verbose If TRUE extra information is printed on the terminal
#' @return A list with the four following elements:
#' \enumerate{
#' \item "FST": a scalar corresponding to the estimate of the genome-wide FST over all the populations
#' \item "snp.FST": a vector containing estimates of SNP-specific FST
#' \item "snp.Q1": a vector containing estimates of the overall within pop. SNP-specific probability of identity
#' \item "snp.Q2": a vector containing estimates of the overall between pop. SNP-specific probability of identity
#' \item "mean.fst" (if nsnp.per.bjack.block>0): genome-wide Fst estimate as the mean over block-jackknife samples (may slight differ from "FST" estimate since it is only computed on SNPs eligible for Block-Jackknife)
#' \item "se.fst" (if nsnp.per.bjack.block>0): standard-error of the genome-wide Fst estimate computed block-jackknife samples
#' \item "fst.bjack.samples" (if nsnp.per.bjack.block>0): a vector containing estimates of the overall between pop. SNP-specific probability of identity
#' \item "sliding.windows.fst" (if sliding.window.size>0): a 4-columns data frame containing information on multi-locus Fst computed for sliding windows of SNPs over the whole genome with i) column with the chromosome/contig of origin of each window; ii) the mid-position of each window; iii) the cumulated mid-position of each window (to facilitate further plotting); and iv) the estimated multi-locus Fst
#' }
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fst=computeFST(pooldata)
#' @export
computeFST<-function(x,method="Anova",nsnp.per.bjack.block=0,sliding.window.size=0,verbose=TRUE){
  if(!(method %in% c("Identity","Anova"))){stop("method should either be Identity or Anova (default)")}
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")} 
  if(method=="Identity"){
    #computing Q1 for each snp/pop
    if(is.countdata(x)){
      snp.Q1=rowMeans((x@refallele.count*(x@refallele.count-1) + (x@total.count-x@refallele.count)*(x@total.count-x@refallele.count-1) )/(x@total.count*(x@total.count-1)),na.rm=T) 
      hat.Q1=mean(snp.Q1,na.rm=T) 
      #computing Q2 for each snp/pop pair
      Q2.countdiff=matrix(0,x@nsnp,x@npops*(x@npops-1)/2)
      Q2.counttot=matrix(0,x@nsnp,x@npops*(x@npops-1)/2)
      cnt=0
      for(i in 1:(x@npops-1)){
        for(j in (i+1):x@npops){
          cnt=cnt+1
          Q2.countdiff[,cnt]=( x@refallele.count[,i]*x@refallele.count[,j] + (x@total.count[,i]-x@refallele.count[,i])*(x@total.count[,j]-x@refallele.count[,j]))
          Q2.counttot[,cnt]=(x@total.count[,i]*x@total.count[,j])
        }
      }
      snp.Q2=rowSums(Q2.countdiff)/rowSums(Q2.counttot)
      rm(Q2.counttot,Q2.countdiff) ; gc()      
      hat.Q2=mean(snp.Q2,na.rm=TRUE)
    }else{#see eq. A39 (Hivert et al., 2018)
      Q1=(x@refallele.readcount*(x@refallele.readcount-1) + (x@readcoverage-x@refallele.readcount)*(x@readcoverage-x@refallele.readcount-1) )/(x@readcoverage*(x@readcoverage-1))
      Q1 = (1/(matrix(1,x@nsnp,x@npools) %*% diag(x@poolsizes-1)))*(Q1 %*% diag(x@poolsizes) - 1) 
      lambdaj=x@poolsizes*(x@poolsizes-1)
      lambdaj=lambdaj/sum(lambdaj)
      snp.Q1=rowSums(Q1%*%diag(lambdaj))
      hat.Q1=mean(snp.Q1,na.rm=T) 
      rm(Q1) ; gc()
      #computing Q2 for each snp/pop pair
      Q2=matrix(0,x@nsnp,x@npools*(x@npools-1)/2)
      omegajj=rep(0,x@npools*(x@npools-1)/2)   
      cnt=0
      for(i in 1:(x@npools-1)){
        for(j in (i+1):x@npools){
          cnt=cnt+1
          omegajj[cnt]=x@poolsizes[i]*x@poolsizes[j]
          Q2[,cnt]=( x@refallele.readcount[,i]*x@refallele.readcount[,j] +
                       (x@readcoverage[,i]-x@refallele.readcount[,i])*(x@readcoverage[,j]-x@refallele.readcount[,j]))/(x@readcoverage[,i]*x@readcoverage[,j])
        }
      }
      snp.Q2=rowSums(Q2%*%diag(omegajj/sum(omegajj)))
      hat.Q2=mean(snp.Q2,na.rm=TRUE)
      rm(Q2) ; gc()
    }
    rslt=list(FST=(hat.Q1-hat.Q2)/(1-hat.Q2),snp.Q1=snp.Q1,snp.Q2=snp.Q2,snp.FST=(snp.Q1-snp.Q2)/(1-snp.Q2) )
    rm(snp.Q1,snp.Q2) ; gc()
  }
  
  if (method=="Anova"){
    if(is.countdata(x)){#eq 5.2 in Weir 1996
      SumNi=rowSums(x@total.count)
      Nic=x@total.count-(x@total.count**2)/SumNi
      Nc=rowSums(Nic)/(x@npops-1)
      MSG=(rowSums(x@refallele.count*(x@total.count-x@refallele.count)/x@total.count)) /(SumNi-x@npops)
      PA=rowSums(x@refallele.count)/SumNi
      MSP=(rowSums(x@total.count*((x@refallele.count/x@total.count-PA)**2)))/(x@npops-1)
      F_ST=(MSP-MSG)/(MSP+(Nc-1)*MSG)
      F_ST_multi=mean(MSP-MSG,na.rm=T)/mean(MSP+(Nc-1)*MSG,na.rm=T)
      rslt <- list(snp.FST = F_ST,snp.Q1 = 1 - MSG*2,snp.Q2 = 1 - MSG*2 - 2*(MSP - MSG) / Nc,FST = F_ST_multi)
      #le facteur 2 dans le calcul de snp.Q1 et snp.Q2 vient du mode de calcul simplifie des MSP et MSG (en bi-allelique) avec les formules a la Weir
      #On peut retomber strictement sur les formules de Rousset 2007 (etendues dans Hivert et al. par Renaud pour le cas PoolSeq) qui matchent parfaitement l'approche poolseq ci-desous en faisant
      # S1 <- rowSums(x@total.count) ; S2 <- rowSums(x@total.count**2)
      # n_c=(S1-S2/S1)/(x@npops-1)
      # YY_alt=x@total.count - x@refallele.count
      # SSI = rowSums(x@refallele.count -  x@refallele.count^2 / x@total.count +  YY_alt - YY_alt^2 / x@total.count,na.rm = TRUE)
      # Ici on calcule les identites intra avec la formule:
      #  Q1=y - y^2/n + (n-y) - (n-y)^2/n
      #    =n - ((y+(n-y))^2 - 2y(n-y))/n
      #    =n - (n^2 - 2y(n-y))/n
      #    =2y(n-y)/n [utilise dans le calcul du MSG ci-dessus]
      # SSP=rowSums(x@total.count * ((x@refallele.count / x@total.count) - (rowSums(x@refallele.count) / S1))^2 + x@total.count * ((YY_alt / x@total.count) - (rowSums(YY_alt) / S1))^2,na.rm = TRUE)
      # MSI=SSI/(S1-x@npops) ; MSP=SSP/(x@npops-1)
      # F_ST=(MSP-MSI)/(MSP+(n_c-1)*MSI)  #cf eq. 28A28 de Rousset (2007) en factorisant numerateur et denominateur par (S1-ns)*(ns-1) pour retrouver MSI et MSP
      # F_ST_multi=mean(MSP-MSI,na.rm=T)/mean(MSP+(n_c-1)*MSI,na.rm=T)
      # snp.Q1 = 1 - MSI  #cf eq 28A21 de Rousset 2007 en sommant sur tous les allele (i.e., SumPi_k=1)
      # snp.Q2 = 1 - MSI - (MSP - MSI) / n_c
      snpNc=Nc #useful for blockjackknife ou sliding wnidows
      rm(F_ST,MSP,MSG) ; gc()
    }else{#Hivert et al. 2018
      mtrx.n_i <- matrix(x@poolsizes,nrow = x@nsnp,ncol = x@npools,byrow = TRUE)
      C_1 <- rowSums(x@readcoverage)
      C_2 <- rowSums(x@readcoverage^2)
      D_2 <- rowSums(x@readcoverage / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i,na.rm = TRUE)
      D_2.star <- rowSums(x@readcoverage * (x@readcoverage / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i),na.rm = TRUE) / C_1
      n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
      rm(C_2,mtrx.n_i) ; gc()
      #    YY_alt <- x@readcoverage - x@refallele.readcount
      #    SSI <- rowSums(x@refallele.readcount - x@refallele.readcount^2 / x@readcoverage + YY_alt - YY_alt^2 / x@readcoverage,na.rm = TRUE)
      SSI <- 2*rowSums(x@refallele.readcount*(x@readcoverage-x@refallele.readcount)/x@readcoverage,na.rm=TRUE)
      SA<-rowSums(x@refallele.readcount)/C_1
      SSP<-2*rowSums(x@readcoverage*((x@refallele.readcount/x@readcoverage-SA)**2),na.rm=TRUE)
      #    SSP <- rowSums(x@readcoverage * ((x@refallele.readcount / x@readcoverage) - (rowSums(x@refallele.readcount) / C_1))^2 + x@readcoverage * ((YY_alt / x@readcoverage) - (rowSums(YY_alt) / C_1))^2,na.rm = TRUE)
      #    rm(YY_alt,mtrx.n_i) ; gc()
      MSI <- SSI / (C_1 - D_2)
      MSP <- SSP / (D_2 - D_2.star)
      rm(C_1,D_2,D_2.star) ; gc()
      F_ST <- (MSP - MSI)  / (MSP + (n_c - 1) * MSI)
      F_ST_multi <- sum(MSP - MSI,na.rm=T)  / sum(MSP + (n_c- 1) * MSI,na.rm=T)
      rslt <- list(snp.FST = F_ST,snp.Q1 = 1 - MSI,snp.Q2 = 1 - MSI - (MSP - MSI) / n_c,FST = F_ST_multi)
      snpNc=n_c #useful for blockjackknife ou sliding wnidow
      rm(F_ST,MSI,MSP,n_c) ; gc()
    }
    snpNc[is.na(rslt$snp.Q1)|is.na(rslt$snp.Q2)]=NA #important pour comptages des SNPs contribuant a Fst dans calcul multilocus (blockjackknife ou sliding window)
    }
  #######
  ##Remarque calcul Fst multilocus (blockjacknife+sliding window) dans le cas Anova à partir des snp.Q1 et snp.Q2
  ##On utilise: Q1-Q2=(MSP-MSI)/nc et (1-Q2)=(MSP+(nc-1)*MSI)/nc
  ##      =>    Fst_multi=sum(nc_i*(Q1_i-Q2_i))/sum(nc_i-nc_i*Q2_i) = sum(MSP_i-MSI_i)/sum(MSP_i+(nc_i-1)*MSP_i) 
  ##Cett procudure tient automatiquement compte des marqueurs NA sur Q2 ou Q1 (sous reserve que snpNC soit a NA)
  ##Dans le cas Identity: il faut aussi compter les NA (sum(Q1-Q2)/(sum(1-Q2)) avec sum(1-Q2)=NmrkNonNA - sum(Q2,na.rm=T)      
  #######
  if(nsnp.per.bjack.block>0){
    if(verbose){cat("Starting Block-Jackknife sampling\n")}
    bjack.blocks=generate.jackknife.blocks(x,nsnp.per.bjack.block,verbose=verbose)
    tmp.idx.sel=!is.na(bjack.blocks$snp.block.id)
    tmp.snp.block.id=bjack.blocks$snp.block.id[!is.na(bjack.blocks$snp.block.id)]
    if(method=="Anova"){
      tmp.sampled.q1=as.vector(by(rslt$snp.Q1[tmp.idx.sel]*snpNc[tmp.idx.sel],tmp.snp.block.id,sum,na.rm=T))     
      tmp.sampled.q1=(sum(tmp.sampled.q1)-tmp.sampled.q1)
      tmp.sampled.q2=as.vector(by(rslt$snp.Q2[tmp.idx.sel]*snpNc[tmp.idx.sel],tmp.snp.block.id,sum,na.rm=T))
      tmp.sampled.q2=(sum(tmp.sampled.q2)-tmp.sampled.q2)      
      tmp.sampled.sumnc=as.vector(by(snpNc[tmp.idx.sel],tmp.snp.block.id,sum,na.rm=T))
      tmp.sampled.sumnc=sum(tmp.sampled.sumnc)-tmp.sampled.sumnc   
      sampled.fst=(tmp.sampled.q1-tmp.sampled.q2) / (tmp.sampled.sumnc-tmp.sampled.q2)  
    }else{
      tmp.sampled.q1=as.vector(by(rslt$snp.Q1[tmp.idx.sel],tmp.snp.block.id,sum,na.rm=T))
      tmp.sampled.q1=(sum(tmp.sampled.q1)-tmp.sampled.q1)
      tmp.sampled.q2=as.vector(by(rslt$snp.Q2[tmp.idx.sel],tmp.snp.block.id,sum,na.rm=T))
      tmp.sampled.q2=(sum(tmp.sampled.q2)-tmp.sampled.q2)
      tmp.snp.cnt=as.vector(by(1-(is.na(rslt$snp.Q1[tmp.idx.sel])|is.na(rslt$snp.Q2[tmp.idx.sel])),tmp.snp.block.id,sum,na.rm=T))
      tmp.snp.cnt=sum(tmp.snp.cnt)-tmp.snp.cnt   
      sampled.fst=(tmp.sampled.q1-tmp.sampled.q2) / (tmp.snp.cnt-tmp.sampled.q2)  
    }    
    rslt[["mean.fst"]]=mean(sampled.fst)
    rslt[["se.fst"]]=sd(sampled.fst)*(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks) #By definition: se(hat(f))={Sum((fi-mean_f)^2)*(n-1)/n}^0.5={var_f*(n-1)*(n-1)/n}^0.5=sd_f*(n-1)/(n^0.5)
    rslt[["fst.bjack.samples"]]=sampled.fst
  }
  
  #######################
  if(sliding.window.size>1){
    if(verbose){cat("Start sliding-window scan\n")}
    det.idx.per.chr=matrix(unlist(by(1:x@nsnp,x@snp.info[,1],range)),ncol=2,byrow=T)
    if(nrow(det.idx.per.chr)==0){#in case all contig names are NA
      cat("Exit function: No chr/contigs available (information on SNP contig name might not have been provided)\n")
    }
    det.idx.per.chr=cbind(det.idx.per.chr,det.idx.per.chr[,2]-det.idx.per.chr[,1]+1) 
    step=floor(sliding.window.size/2)
    all.pos=all.fst=all.chr=all.cumpos=win.size=c()
    tmp.cum=0
    if(verbose){
      n.chr.eval=sum(det.idx.per.chr[,3]>sliding.window.size)
      cat(n.chr.eval,"chromosomes scanned (with more than",sliding.window.size,"SNPs)\n")
      pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.chr.eval, clear = FALSE, width= 60)
    #  pb <- txtProgressBar(min = 0, max = n.chr.eval, style = 3)
      tmp.cnt=0
    }
    all.snp.pos=x@snp.info[,2]
    for(i in 1:nrow(det.idx.per.chr)){
      if(det.idx.per.chr[i,3]>sliding.window.size){
        tmp.sel=det.idx.per.chr[i,1]:det.idx.per.chr[i,2]
        tmp.pos=floor(rollmean(all.snp.pos[tmp.sel],k=sliding.window.size))
        retained.pos=seq(1,length(tmp.pos),step)
        tmp.pos=tmp.pos[retained.pos]
        all.pos=c(all.pos,tmp.pos)
        all.cumpos=c(all.cumpos,tmp.pos+tmp.cum)
        tmp.cum=max(all.cumpos)
        all.chr=c(all.chr,rep(x@snp.info[tmp.sel[1],1],length(tmp.pos)))
        #calcul fst
        if(method=="Anova"){
          qq1=rollsum(rslt$snp.Q1[tmp.sel]*snpNc[tmp.sel],k=sliding.window.size,na.rm=T)[retained.pos]
          qq2=rollsum(rslt$snp.Q2[tmp.sel]*snpNc[tmp.sel],k=sliding.window.size,na.rm=T)[retained.pos]  
          tmp.sumnc=rollsum(snpNc[tmp.sel],k=sliding.window.size,na.rm=T)[retained.pos]
          tmp.fst=(qq1-qq2) / (tmp.sumnc-qq2)   
        }else{
          qq1=rollsum(rslt$snp.Q1[tmp.sel],k=sliding.window.size,na.rm=T)[retained.pos]
          qq2=rollsum(rslt$snp.Q2[tmp.sel],k=sliding.window.size,na.rm=T)[retained.pos]  
          tmp.snp.cnt=rollsum(1-(is.na(rslt$snp.Q1[tmp.sel])|is.na(rslt$snp.Q2[tmp.sel])),k=sliding.window.size)[retained.pos]  
          tmp.fst=(qq1-qq2) / (tmp.snp.cnt-qq2)  
        }    
        all.fst=c(all.fst,tmp.fst)
        win.size=c(win.size,all.snp.pos[tmp.sel[retained.pos]+sliding.window.size-1]-all.snp.pos[tmp.sel[retained.pos]])
        if(verbose){pb$tick()} #tmp.cnt=tmp.cnt+1 ; setTxtProgressBar(pb, tmp.cnt)
      }
    }
    if(verbose){
      cat("\nAverage (min-max) Window Sizes",round(mean(win.size*1e-3),1),"(",round(min(win.size*1e-3),1),"-",round(max(win.size*1e-3),1),") kb\n")
      pb$terminate()}#close(pb)}
    rslt[["sliding.windows.fst"]]=data.frame(Chr=all.chr,Position=all.pos,CumulatedPosition=all.cumpos,MultiLocusFst=all.fst,stringsAsFactors=FALSE)
  }
  
  return(rslt)
}
