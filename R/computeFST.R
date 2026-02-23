#' Compute Fst from Pool-Seq data or Count data
#' @param x A pooldata object containing Pool-Seq information or countdata object containing allele counts information
#' @param method Either "Anova" (default method) or "Identity" (relying on unbiased estimators of Probability of Identity within and across pairs of pools/populations)
#' @param struct Vector of length equal to the number of pop. sample that give the pop. sample group name of index (i.e., structure) 
#' @param weightpid When method="Identity", if TRUE weighting averages of pop. Q1 and pairwise Q2 are performed (see eq. A46 and A47 in Hivert et al., 2018 for PoolSeq and Rousset 2007 for count data) to compute overall Q1 and Q2. If not, unweighted averages are performed. 
#' @param nsnp.per.bjack.block Number of consecutive SNPs within a block for block-jackknife (default=0, i.e., no block-jackknife sampling) 
#' @param sliding.window.size Number of consecutive SNPs within a window for multi-locus computation of Fst over sliding window with half-window size step (default=0, i.e., no sliding-window scan) 
#' @param verbose If TRUE extra information is printed on the terminal
#' @return A list with the four following elements:
#' \enumerate{
#' \item "FST": estimate of genome-wide Fst over all the populations. The element is a vector with 5 elements corresponding to i) the estimated value over all SNPs; ii) the block-jackknife mean; iii)  the block-jackknife s.e.; iv) the lower; and v) the upper bound of the 95% (block-jackknife) Confidence Interval estimates. If nsnp.per.bjack.block=0, only the estimated value is given (other elements are set to NA).
#' \item "FSG": under the hierarchical Fst model (i.e., when struct vector is non-null); estimates estimate of genome-wide within-group differentiation (Fsg). The element is a vector with 5 elements corresponding to i) the estimated value over all SNPs; ii) the block-jackknife mean; iii)  the block-jackknife s.e.; iv) the lower; and v) the upper bound of the 95% (block-jackknife) Confidence Interval estimates. If nsnp.per.bjack.block=0, only the estimated value is given (other elements are set to NA).
#' \item "FGT": under the hierarchical Fst model (i.e., when struct vector is non-null); estimates estimate of genome-wide between-group differentiation (Fgt). The element is a vector with 5 elements corresponding to i) the estimated value over all SNPs; ii) the block-jackknife mean; iii)  the block-jackknife s.e.; iv) the lower; and v) the upper bound of the 95% (block-jackknife) Confidence Interval estimates. If nsnp.per.bjack.block=0, only the estimated value is given (other elements are set to NA).
#' \item "snp.Fstats": a data frame containing SNP-specific estimates of Fst and also under the hierarchical (i.e., when struct vector is non-null) SNP-specific estimates Fsg and Fgt
#' \item "snp.Q": a data frame containing SNP-specific estimates of Q1 (within-population) and Q2 (between-population) probability of identity and also under the hierarchical (i.e., when struct vector is non-null) SNP-specific estimates of Q3, the probability of identity between populations from different groups (under this model Q2 is then the Pid between populations from the same group).
#' \item "sliding.windows.fvalues" (if sliding.window.size>0): a 4 or 6 (under hierarchical Fst model) column data frame containing information on multi-locus Fst (and Fsg and Fgt under the hierarchical Fst model) computed for sliding windows of SNPs over the whole genome with i) column with the chromosome/contig of origin of each window; ii) the mid-position of each window; iii) the cumulated mid-position of each window (to facilitate further plotting); iv) the estimated multi-locus Fst; and under the hierarchical Fst model v) the estimated multi-locus Fsg and ; vi) the estimated multi-locus Fgt
#' }
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fst=computeFST(pooldata)
#'  res.hierfst=computeFST(pooldata,struct=c(rep("A",5),rep("B",7),rep("C",3)))
#' @export
computeFST <- function(x,method = "Anova",struct = NULL,weightpid=FALSE,nsnp.per.bjack.block=0,sliding.window.size=0,verbose=TRUE) {
  #check input objects (and retrieve nbr of pop samples)
  if (!(method %in% c("Identity","Anova"))) {stop("method should either be Identity or Anova (default)")}
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")} 
  if(is.countdata(x)){nbr.pops <- x@npops}else{nbr.pops <- x@npools}
  if(nbr.pops<2){stop("At least 2 pop. samples required to compute Fst\n")} #may arrive after pooldata.subset with a single pop
  if(is.null(struct)){
    compute.hierfstat=FALSE    
  }else{
    #check structure vector  
    if(length(struct)!=nbr.pops){
      stop("Misspecified structure vector (struct argument): must the same length as the number of pop. samples in the countdata or pooldata object.\n")}
    nbr.grps <- length(sort(unique(struct)))
    if(nbr.grps<2){stop("A single group was found for all sampled pop in the structure vector. At least 2 groups must be declared to compute hierarchical Fst!\n")}
    grp.code=1:nbr.grps ; names(grp.code)=unique(struct)
    struct.code=as.numeric(grp.code[struct])
    compute.hierfstat=TRUE 
  }
  
  if(method=="Anova" & weightpid){cat("NOTE: weightpid is not relevant when method is Anova\n")}
  
  if(compute.hierfstat){
    if(verbose){cat("Computation of hierarchical Fst (Gautier et al., 2024)\n")}
    if(verbose){
      cat(nbr.grps,"groups of pop samples declared in struct object:\n")
      print(grp.code)
    }
  }else{
    if(verbose){cat("Computation of Fst (Hivert et al., 2018)\n")}    
  }
  ########
  ##prepare output object
  rslt=list()
  tmp.vect=rep(NA,5) ; names(tmp.vect)=c("Estimate","bjack mean","bjack s.e.","CI95inf","CI95sup")
  rslt$Fst=tmp.vect
  if(compute.hierfstat){
    rslt$snp.Q=rslt$snp.Fstats=as.data.frame(matrix(NA,x@nsnp,3))
    colnames(rslt$snp.Q)=c("Q1","Q2","Q3")
    colnames(rslt$snp.Fstats)=c("Fsg","Fgt","Fst") 
    rslt$Fsg=rslt$Fgt=tmp.vect
  }else{
    rslt$snp.Q=as.data.frame(matrix(NA,x@nsnp,2))
    rslt$snp.Fstats=as.data.frame(matrix(NA,x@nsnp,1))
    colnames(rslt$snp.Q)=c("Q1","Q2")   
    colnames(rslt$snp.Fstats)=c("Fst")  
  }
  ########################
  if(method=="Identity"){
    if(verbose){cat("Computing SNP-specific Q1 (Identity estimator)\n")}
    if(is.countdata(x)){
      if(weightpid){
       rslt$snp.Q[,1]=.compute_snpQ1rw(x@refallele.count,x@total.count,rep(1,nbr.pops),rep(1,nbr.pops),FALSE,verbose=verbose)
      }else{
       rslt$snp.Q[,1]=.compute_snpQ1(x@refallele.count,x@total.count,rep(1,nbr.pops),verbose=verbose) 
      }
    }else{
      if(weightpid){
        rslt$snp.Q[,1]=.compute_snpQ1rw(x@refallele.readcount,x@readcoverage,x@poolsizes/(x@poolsizes-1),x@poolsizes,TRUE,verbose=verbose)
      }else{
        rslt$snp.Q[,1]=.compute_snpQ1(x@refallele.readcount,x@readcoverage,x@poolsizes/(x@poolsizes-1),verbose=verbose) 
      }
    }
    pairs.idx=t(combn(nbr.pops,2))
    if(compute.hierfstat){
      pairs.within=which(struct.code[pairs.idx[,1]]==struct.code[pairs.idx[,2]])
      #computing pairwise within group 
      if(verbose){cat("Computing SNP-specific within-group Q2 (Identity estimator)\n")}
      if(is.countdata(x)){
        if(weightpid){
          rslt$snp.Q[,2]=.compute_snpQ2rw(x@refallele.count,x@total.count,pairs.idx[pairs.within,]-1,rep(1,nbr.pops),FALSE,verbose=verbose) 
        }else{
          rslt$snp.Q[,2]=.compute_snpQ2(x@refallele.count,x@total.count,pairs.idx[pairs.within,]-1,verbose=verbose)
        }
      }else{
        if(weightpid){
          rslt$snp.Q[,2]=.compute_snpQ2rw(x@refallele.readcount,x@readcoverage,pairs.idx[pairs.within,]-1,x@poolsizes,TRUE,verbose=verbose) 
        }else{
         rslt$snp.Q[,2]=.compute_snpQ2(x@refallele.readcount,x@readcoverage,pairs.idx[pairs.within,]-1,verbose=verbose)
        }
      }
      #computing pairwise across group
      if(verbose){cat("Computing SNP-specific Q3 i.e. between-group Q2 (Identity estimator)\n")}
      if(is.countdata(x)){
        if(weightpid){
          rslt$snp.Q[,3]=.compute_snpQ2rw(x@refallele.count,x@total.count,pairs.idx[-pairs.within,]-1,rep(1,nbr.pops),FALSE,verbose=verbose) 
        }else{
          rslt$snp.Q[,3]=.compute_snpQ2(x@refallele.count,x@total.count,pairs.idx[-pairs.within,]-1,verbose=verbose)
        }        
       }else{
         if(weightpid){
           rslt$snp.Q[,3]=.compute_snpQ2rw(x@refallele.readcount,x@readcoverage,pairs.idx[-pairs.within,]-1,x@poolsizes,TRUE,verbose=verbose) 
         }else{
           rslt$snp.Q[,3]=.compute_snpQ2(x@refallele.readcount,x@readcoverage,pairs.idx[-pairs.within,]-1,verbose=verbose)
         }
       }
      num.fsg=(rslt$snp.Q$Q1 - rslt$snp.Q$Q2)
      den.fsg=1. - rslt$snp.Q$Q2
      num.fgt=rslt$snp.Q$Q2 - rslt$snp.Q$Q3
      num.fst=rslt$snp.Q$Q1 - rslt$snp.Q$Q3
      den.fst=1. - rslt$snp.Q$Q3 #same as den.fgt	    
    }else{
      #computing pairwise within group 
      if(verbose){cat("Computing SNP-specific within-group Q2 (Identity estimator)\n")}
      if(is.countdata(x)){
        if(weightpid){
          rslt$snp.Q[,2]=.compute_snpQ2rw(x@refallele.count,x@total.count,pairs.idx-1,rep(1,nbr.pops),FALSE,verbose=verbose) 
        }else{
          rslt$snp.Q[,2]=.compute_snpQ2(x@refallele.count,x@total.count,pairs.idx-1,verbose=verbose)
        }
      }else{
        if(weightpid){
          rslt$snp.Q[,2]=.compute_snpQ2rw(x@refallele.readcount,x@readcoverage,pairs.idx-1,x@poolsizes,TRUE,verbose=verbose) 
        }else{
          rslt$snp.Q[,2]=.compute_snpQ2(x@refallele.readcount,x@readcoverage,pairs.idx-1,verbose=verbose)
        }        
      }
      num.fst=rslt$snp.Q$Q1 - rslt$snp.Q$Q2
      den.fst=1. - rslt$snp.Q$Q2 #same as den.fgt
    }
  }
  ########################
  if(method=="Anova"){
    if(compute.hierfstat){
      if(verbose){cat("Computing SNP-specific Q1, Q2 and Q3 (Anova estimator)\n")}
      if(is.countdata(x)){
        tmp.aov=.compute_snpHierFstAov(x@refallele.count,x@total.count,rep(-1,nbr.pops),struct.code-1,verbose=verbose)
      }else{
        tmp.aov=.compute_snpHierFstAov(x@refallele.readcount,x@readcoverage,x@poolsizes,struct.code-1,verbose=verbose)
      }
      #tmp.aov: MSI, MSP, MSG, nc, nc_p, nc_pp
      rslt$snp.Q[,1]=1.-tmp.aov[,1]                                         #1-MSI
      rslt$snp.Q[,2]=rslt$snp.Q[,1] - (tmp.aov[,2]-tmp.aov[,1])/tmp.aov[,4] #1-MSI-(MSP - MSI) / n_c
      rslt$snp.Q[,3]=rslt$snp.Q[,2] - 
        ((tmp.aov[,3]-tmp.aov[,1])+(tmp.aov[,2]-tmp.aov[,1])*tmp.aov[,6]/tmp.aov[,4])/tmp.aov[,5]           
      #Q3=1.0 - MSI - (MSP - MSI) / nc - (MSG - MSI + (MSP - MSI) * nc_pp / nc) / nc_p
      #important de ne pas passer par les Q1, Q2 et Q3 pour garder le facteur nc le calcul multilocus
      num.fsg=(tmp.aov[,2]-tmp.aov[,1])      #=(MSP - MSI)
      den.fsg=(tmp.aov[,2]+(tmp.aov[,4]-1.)*tmp.aov[,1])   #=(MSP + (nc - 1.0) * MSI)
      num.fgt=tmp.aov[,4]*(tmp.aov[,3]-tmp.aov[,1]) + tmp.aov[,6]*(tmp.aov[,2]-tmp.aov[,1])
      #num.fgt=(nc * (MSG - MSI) + nc_pp * (MSP - MSI))
      num.fst=tmp.aov[,4]*(tmp.aov[,3]-tmp.aov[,1]) + (tmp.aov[,5]+tmp.aov[,6])*(tmp.aov[,2]-tmp.aov[,1])
      #	num.fst=(nc * (MSG - MSI) + (nc_p + nc_pp) * (MSP - MSI))
      den.fst=num.fst + tmp.aov[,4]*tmp.aov[,5]*tmp.aov[,1]
      #den.fgt=den.fst=(nc * (MSG - MSI) + (nc_p + nc_pp) * (MSP - MSI) + nc * nc_p * MSI)
    }else{
      if(verbose){cat("Computing SNP-specific Q1 and Q2 (Anova estimator)\n")}
      if(is.countdata(x)){
        tmp.aov=.compute_snpFstAov(x@refallele.count,x@total.count,rep(-1,nbr.pops),verbose=verbose)
      }else{
        tmp.aov=.compute_snpFstAov(x@refallele.readcount,x@readcoverage,x@poolsizes,verbose=verbose)
      }
      rslt$snp.Q[,1]=1.-tmp.aov[,1]                                         #1-MSG
      rslt$snp.Q[,2]=rslt$snp.Q[,1] - (tmp.aov[,2]-tmp.aov[,1])/tmp.aov[,3] #1 - MSG - (MSP - MSG) / n_c       
      num.fst=tmp.aov[,2]-tmp.aov[,1]                          #MSP-MSG :
      den.fst=tmp.aov[,2]+(tmp.aov[,3]-1)*tmp.aov[,1]          #MSP + (n_c- 1) * MSG 
      #important de ne pas passer par les Q1 et Q2 pour garder le facteur nc le calcul multilocus
    }
    rm(tmp.aov) ; gc()
  }
  ########################
  #multi-locus estimates
  keep=rowSums(is.na(rslt$snp.Q))==0
  if(compute.hierfstat){
    rslt$snp.Fstats$Fsg=num.fsg / den.fsg  
    rslt$snp.Fstats$Fgt=num.fgt / den.fst
    rslt$snp.Fstats$Fst=num.fst / den.fst
    rslt$Fsg[1]=sum(num.fsg[keep])/sum(den.fsg[keep]) 
    rslt$Fgt[1]=sum(num.fgt[keep])/sum(den.fst[keep]) 
    rslt$Fst[1]=sum(num.fst[keep])/sum(den.fst[keep]) 
  }else{
    rslt$snp.Fstats$Fst=num.fst / den.fst
    rslt$Fst[1]=sum(num.fst[keep])/sum(den.fst[keep])   
  }
  
  #######
  #block jackknife      
  #######
  if(nsnp.per.bjack.block>0){
    if(verbose){cat("Starting Block-Jackknife sampling\n")}
    bjack.blocks=generate.jackknife.blocks(x,nsnp.per.bjack.block,verbose=verbose)
    keep=keep & !is.na(bjack.blocks$snp.block.id)
    tmp.snp.block.id=bjack.blocks$snp.block.id[keep]-1 #0 indexed for cpp
    
    num.fst.blk.contrib=.block_sum(num.fst[keep],tmp.snp.block.id)
    den.fst.blk.contrib=.block_sum(den.fst[keep],tmp.snp.block.id)
    den.fst.blk.contrib=den.fst.blk.contrib-sum(den.fst.blk.contrib) #keep it for fgt (if hierfst)     
    sampled.fst=(num.fst.blk.contrib-sum(num.fst.blk.contrib))/den.fst.blk.contrib
    rm(num.fst.blk.contrib) ; gc()    #keep   den.fst.blk.contrib for fgt if hierfstat
    tmp.fact=(bjack.blocks$nblocks-1)/sqrt(bjack.blocks$nblocks)
    rslt$Fst[2:3]=c(mean(sampled.fst),sd(sampled.fst)*tmp.fact)  
    rm(sampled.fst) ; gc()
    rslt$Fst[4:5]=rslt$Fst[2]+c(-1.96,1.96)*rslt$Fst[3]   
    if(compute.hierfstat){
      num.fsg.blk.contrib=.block_sum(num.fsg[keep],tmp.snp.block.id)
      den.fsg.blk.contrib=.block_sum(den.fsg[keep],tmp.snp.block.id)
      sampled.fsg=(num.fsg.blk.contrib-sum(num.fsg.blk.contrib))/(den.fsg.blk.contrib-sum(den.fsg.blk.contrib))
      rm(num.fsg.blk.contrib,den.fsg.blk.contrib) ; gc()
      rslt$Fsg[2:3]=c(mean(sampled.fsg),sd(sampled.fsg)*tmp.fact)
      rm(sampled.fsg) ; gc()
      rslt$Fsg[4:5]=rslt$Fsg[2]+c(-1.96,1.96)*rslt$Fsg[3]   
      
      num.fgt.blk.contrib=.block_sum(num.fgt[keep],tmp.snp.block.id)
      sampled.fgt=(num.fgt.blk.contrib-sum(num.fgt.blk.contrib))/den.fst.blk.contrib
      rm(num.fgt.blk.contrib) ; gc()
      rslt$Fgt[2:3]=c(mean(sampled.fgt),sd(sampled.fgt)*tmp.fact)
      rm(sampled.fgt) ; gc()      
      rslt$Fgt[4:5]=rslt$Fgt[2]+c(-1.96,1.96)*rslt$Fgt[3]
    }
    rm(den.fst.blk.contrib) ; gc()
  }
  #######
  #sliding windows     
  #######
  if(sliding.window.size>1){
    keep=rowSums(is.na(rslt$snp.Q))==0 #may have been modified by blockjackknife
    if(sum(keep)<1){stop("Error no SNP eligible\n")}
    if(verbose){cat("Start sliding-window scan\n")}
    det.idx.per.chr=matrix(unlist(by(1:sum(keep),x@snp.info[keep,1],range)),ncol=2,byrow=T)
    if(nrow(det.idx.per.chr)==0){#in case all contig names are NA
      cat("Exit function: No chr/contigs available (information on SNP contig name might not have been provided)\n")
    }
    det.idx.per.chr=cbind(det.idx.per.chr,det.idx.per.chr[,2]-det.idx.per.chr[,1]+1) 
    step=floor(sliding.window.size/2)
    all.pos=all.startpos=all.endpos=all.fsg=all.fgt=all.fst=all.chr=all.cumpos=win.size=c()
    tmp.cum=0
    if(verbose){
      n.chr.eval=sum(det.idx.per.chr[,3]>sliding.window.size)
      cat(n.chr.eval,"chromosomes scanned (with more than",sliding.window.size,"SNPs)\n")
      pb <- progress_bar$new(format = " [:bar] :percent eta: :eta",total = n.chr.eval, clear = FALSE, width= 60)
      tmp.cnt=0
    }
    all.snp.pos=x@snp.info[keep,2]
    for(i in 1:nrow(det.idx.per.chr)){
      if(det.idx.per.chr[i,3]>sliding.window.size){
        tmp.sel=det.idx.per.chr[i,1]:det.idx.per.chr[i,2] #idx all SNPs for ith chrom.
        tmp.win.start=seq(1,det.idx.per.chr[i,3]-sliding.window.size,step)
        tmp.win.end=tmp.win.start+sliding.window.size-1
        tmp.sel=tmp.sel[1:rev(tmp.win.end)[1]] #remove remaining SNPs after the last window
        tmp.nwin=length(tmp.win.start)
        tmp.snp.win.idx=cbind(tmp.win.start,tmp.win.end)-1 #0-indexed (for cpp)
        tmp.den.fst=.block_sum2(den.fst[tmp.sel],tmp.snp.win.idx) #to keep for fgt if needed
        all.fst=c(all.fst,.block_sum2(num.fst[tmp.sel],tmp.snp.win.idx)/tmp.den.fst) 
        if(compute.hierfstat){
          all.fsg=c(all.fsg,.block_sum2(num.fsg[tmp.sel],tmp.snp.win.idx)/.block_sum2(den.fsg[tmp.sel],tmp.snp.win.idx))  
          all.fgt=c(all.fgt,.block_sum2(num.fgt[tmp.sel],tmp.snp.win.idx)/tmp.den.fst)
        }
        all.startpos=c(all.startpos,all.snp.pos[tmp.sel][tmp.win.start])
        all.endpos=c(all.endpos,all.snp.pos[tmp.sel][tmp.win.end])
        tmp.pos=(all.snp.pos[tmp.sel][tmp.win.start]+all.snp.pos[tmp.sel][tmp.win.end])/2
        all.pos=c(all.pos,tmp.pos)
        all.cumpos=c(all.cumpos,tmp.pos+tmp.cum)
        tmp.cum=max(all.cumpos)
        all.chr=c(all.chr,rep(x@snp.info[keep,1][tmp.sel[1]],tmp.nwin))       
        tmp.winsize=all.snp.pos[tmp.sel][tmp.win.end]-all.snp.pos[tmp.sel][tmp.win.start]
        win.size=c(win.size,tmp.winsize)
        if(verbose){pb$tick()} #tmp.cnt=tmp.cnt+1 ; setTxtProgressBar(pb, tmp.cnt)
      }
    }
    
    if(verbose){
      cat("\nAverage (min-max) Window Sizes",round(mean(win.size*1e-3),1),"(",round(min(win.size*1e-3),1),"-",round(max(win.size*1e-3),1),") kb\n")
      pb$terminate()}#close(pb)}
    if(compute.hierfstat){
      rslt$sliding.windows.fvalues=data.frame(Chr=all.chr,Start=all.startpos,End=all.endpos,MidPos=all.pos,CumMidPos=all.cumpos,
                                              MultiLocusFsg=all.fsg,MultiLocusFgt=all.fgt,MultiLocusFst=all.fst,stringsAsFactors=FALSE)
    }else{
      rslt$sliding.windows.fvalues=data.frame(Chr=all.chr,Start=all.startpos,End=all.endpos,MidPos=all.pos,CumMidPos=all.cumpos,
                                              MultiLocusFst=all.fst,stringsAsFactors=FALSE)        
    }
  }
  return(rslt)	
}
