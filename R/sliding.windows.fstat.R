#' Compute sliding window estimates of F-statistics or ratio of F-statistics over the genome
#' @param x A pooldata object containing Pool-Seq information or a countdata object containing allele count information
#' @param num.pop.idx A vector of length 1 to 4 (depending on num.stat) giving the index (or names) of the populations samples. If num.stat="het", num.pop.idx must be of length 1: num.pop.idx=i specifies the ith pop in x. If num.stat="div", "F2" or "Fst", num.pop.idx must be of length 2: num.pop.idx=c(i,j) specifies the pairs of populations with indexes i and j in x. If num.stat="F3" or "F3star", num.pop.idx must be of length 3 (num.pop.idx=c(i,j,k) specifies the F3(pop_i;pop_j,pop_k) populations triplet). Finally, if num.stat="F4" or "Dstat", num.pop.idx must be of length 4: num.pop.idx=c(i,j,k) specifies the F4(pop_i,pop_j;pop_k,pop_l) populations quadruplet i.e. the computed (numerator) statistic computed is (F2(pop_i,pop_k)-F2(pop_i,pop_l)-F2(pop_j,pop_k)+F2(pop_j,pop_l))/2.
#' @param den.pop.idx A vector of length 1 to 4 (see num.pop.idx description) giving the index of the populations specifying the F-statistic. If NULL, the computed statistic is the one specified by num.pop.idx.
#' @param num.stat the name of the (numerator) stat which must be "het" (1-Q1), "div" (1-Q2), "F2", "Fst", "F3", "F3star", "F4", "Dstat", "Fh", "Fd", or "FdM"
#' @param den.stat the name of the (numerator) stat which must be "het" (1-Q1), "div" (1-Q2), "F2", "Fst", "F3", "F3star", "F4", "Dstat", "Fh", "Fd", or "FdM"
#' @param window.def Either "nsnp" or "bp" to define windows by either a number of SNPs or a size in bp, respectively 
#' @param sliding.window.size A numeric value giving the number of SNPs or the size (in bp) of the windows depending window.def
#' @param window.overlap.fact A numeric value (between 0 and 1) giving the percentage of overlap between consecutive windows (default=0.5)
#' @param bp.start.first.snp When window.def="bp", if TRUE (default) the windowing start at the first SNP position, if FALSE the windowing start at position 1
#' @param verbose If TRUE extra information is printed on the terminal
#' @details Compute sliding window estimates of F-statistics or ratio of F-statistics over the genome.
#' @return A data frame with 7 columns with for each window in a row their i) chromosome/contig of origin; ii) start and iii) end position; iv) the mid-position of each window; v) the cumulated mid-position of each window (to facilitate further plotting); vi) the number of SNPs included in the computation of window value; and vii) the estimated value of the statistic
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}},\code{\link{genobaypass2pooldata}} or \code{\link{genoselestim2pooldata}}. To generate coundata object, see \code{\link{genobaypass2countdata}} or \code{\link{genotreemix2countdata}}.
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#' @export
sliding.windows.fstat<-function(x,num.pop.idx=NULL,den.pop.idx=NULL,num.stat=NULL,den.stat=NULL,window.def=c("nsnp","bp")[1],sliding.window.size=NULL,window.overlap.fact=0.5,bp.start.first.snp=TRUE,verbose=TRUE){
  
  #check input objects (and retrieve nbr of pop samples)
  if(!(is.pooldata(x)) & !(is.countdata(x))){
    stop("Input data are not formatted as valid pooldata (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions) or countdata (see the genobaypass2countdata and genotreemix2countdata) object\n")} 
  if(is.countdata(x)){
    n.pops=x@npops ; pop.names=x@popnames
  }else{
    n.pops=x@npools ; pop.names=x@poolnames
  }
  npop.stats=c(1,2,2,2,3,3,4,4,4,4,4)
  names(npop.stats)=c("het","div","F2","Fst","F3","F3star","F4","Dstat","Fh","Fd","FdM")
  if(is.null(num.pop.idx)){stop("No pop. indexes (or names) given for the numerator stat. (num.pop.idx argument)\n")}
  if(is.null(num.stat)){stop("No numerator f-stat name given (num.stat argument)\n")}  
  if(!(num.stat %in% names(npop.stats))){
    stop(paste0("num.stat must be one of the following names: ",paste(names(npop.stats),collapse=", ")))
  }
  if(length(num.pop.idx)!=npop.stats[num.stat]){
    stop(paste("num.pop.idx argument for",num.stat,"must be of length",npop.stats[num.stat])) 
  }
  #check whether consistent indexes or names for num.pop.idx
  if(is.numeric(num.pop.idx)){
    if(sum(!(num.pop.idx%in%1:n.pops))>0){
      stop("Check num.pop.idx indexes (must be >=1 and <=n.pops included in the countdata or pooldata input object)\n")
    }
  }else{
    if(sum(!(num.pop.idx%in%pop.names))>0){
      stop("Check sample names (detected in num.pop.idx): all samples must be included in the countdata or pooldata input object\nYou may also provide indexes instead\n")
    }
    num.pop.idx=sapply(num.pop.idx, function(i) which(pop.names == i))
  }
  
  if(is.null(den.pop.idx) | is.null(den.stat)){
    denom=FALSE
  }else{
    if(is.null(den.pop.idx)){stop("No pop. indexes given for the denominator stat. (den.pop.idx argument)\n")}
    if(is.null(den.stat)){stop("No denominator f-stat name given (den.stat argument)\n")}  
    if(!(den.stat %in% names(npop.stats))){
      stop(paste0("den.stat must be one of the following names: ",paste(names(npop.stats),collapse=", ")))
    }
    if(length(den.pop.idx)!=npop.stats[den.stat]){
      stop(paste("den.pop.idx argument for",den.stat,"must be of length",npop.stats[den.stat])) 
    }  
    #check whether consistent indexes or names for den.pop.idx
    if(is.numeric(den.pop.idx)){
      if(sum(!(den.pop.idx%in%1:n.pops))>0){
        stop("Check den.pop.idx indexes (must be >=1 and <=n.pops included in the countdata or pooldata input object)\n")
      }
    }else{
      if(sum(!(den.pop.idx%in%pop.names))>0){
        stop("Check sample names (detected in den.pop.idx): all samples must be included in the countdata or pooldata input object\nYou may also provide indexes instead\n")
      }
      den.pop.idx=sapply(den.pop.idx, function(i) which(pop.names == i))      
    }
    denom=TRUE
  }    
  
  if (!(window.def %in% c("nsnp","bp"))){
    stop("window.def should either be nsnp (windows are defined by a number of SNPs, default) or bp (windows are defined by a size in bp)")
  }else{
    winbp=ifelse(window.def=="nsnp",FALSE,TRUE)
  }
  if(is.null(sliding.window.size)){stop("No sliding.window.size provided\n")}
  if(window.overlap.fact<=0 | window.overlap.fact>=1){stop("window.overlap.fact must be between 0 and 1\n")} 
  step=floor(sliding.window.size*(1-window.overlap.fact))
  ###############
  ########internal function to manage computation of the stats (all checks are assumed to be OK)
  ########NOTE: only compute numerator of scaled stats
  ########
  compute_snp_stat_ <- function(statname,pop.idx){
    if(statname=="het"){
      if(is.countdata(x)){
        tmp.w=c(1)
        tmp.refcount=matrix(x@refallele.count[,pop.idx],ncol=1)
        tmp.totcount=matrix(x@total.count[,pop.idx],ncol=1) 
      }else{
        tmp.w=x@poolsizes[pop.idx]/(x@poolsizes[pop.idx]-1)
        tmp.refcount=matrix(x@refallele.readcount[,pop.idx],ncol=1)
        tmp.totcount=matrix(x@readcoverage[,pop.idx],ncol=1)          
      }
      snp.val=1-.compute_snpQ1onepop(tmp.refcount,tmp.totcount,tmp.w)
    } 
    
    if(statname=="div"){
      if(is.countdata(x)){
        snp.val=1-.compute_snpQ2onepair(x@refallele.count[,pop.idx[1]],x@refallele.count[,pop.idx[2]],
                                        x@total.count[,pop.idx[1]],x@total.count[,pop.idx[2]])
      }else{
        snp.val=1-.compute_snpQ2onepair(x@refallele.readcount[,pop.idx[1]],x@refallele.readcount[,pop.idx[2]],
                                        x@readcoverage[,pop.idx[1]],x@readcoverage[,pop.idx[2]])
      }    
    }    
    
    if(statname%in%c("F2","Fst")){
      snp.val=matrix(0,x@nsnp,3) #Q1a Q1b and Q2ab
      for(k in 1:2){
        if(is.countdata(x)){
          tmp.w=c(1)
          tmp.refcount=matrix(x@refallele.count[,pop.idx[k]],ncol=1)
          tmp.totcount=matrix(x@total.count[,pop.idx[k]],ncol=1) 
        }else{
          tmp.w=x@poolsizes[pop.idx[k]]/(x@poolsizes[pop.idx[k]]-1)
          tmp.refcount=matrix(x@refallele.readcount[,pop.idx[k]],ncol=1)
          tmp.totcount=matrix(x@readcoverage[,pop.idx[k]],ncol=1)          
        }          
        snp.val[,k]=.compute_snpQ1onepop(tmp.refcount,tmp.totcount,tmp.w)
      }
      if(is.countdata(x)){
        tmp.refcount=x@refallele.count[,pop.idx];tmp.totcount=x@total.count[,pop.idx]
      }else{
        tmp.refcount=x@refallele.readcount[,pop.idx];tmp.totcount=x@readcoverage[,pop.idx]        
      }           
      snp.val[,3]=.compute_snpQ2onepair(tmp.refcount[,1],tmp.refcount[,2],tmp.totcount[,1],tmp.totcount[,2])
      snp.val=0.5*rowSums(snp.val[,1:2]) - snp.val[,3]
    }
    
    if(statname%in%c("F3","F3star")){
      snp.val=matrix(0,x@nsnp,4) #Q1a Q2bc Q2ab Q2ac
      if(is.countdata(x)){
        tmp.w=c(1)
        tmp.refcount=matrix(x@refallele.count[,pop.idx[1]],ncol=1)
        tmp.totcount=matrix(x@total.count[,pop.idx[1]],ncol=1) 
      }else{
        tmp.w=x@poolsizes[pop.idx[1]]/(x@poolsizes[pop.idx[1]]-1)
        tmp.refcount=matrix(x@refallele.readcount[,pop.idx[1]],ncol=1)
        tmp.totcount=matrix(x@readcoverage[,pop.idx[1]],ncol=1)          
      }      
      snp.val[,1]=.compute_snpQ1onepop(tmp.refcount,tmp.totcount,tmp.w)
      if(is.countdata(x)){
        snp.val[,2]=.compute_snpQ2onepair(x@refallele.count[,pop.idx[2]],x@refallele.count[,pop.idx[3]],
                                          x@total.count[,pop.idx[2]],x@total.count[,pop.idx[3]])
        for(k in 2:3){#ab et ac
          snp.val[,k+1]=.compute_snpQ2onepair(x@refallele.count[,pop.idx[1]],x@refallele.count[,pop.idx[k]],
                                              x@total.count[,pop.idx[1]],x@total.count[,pop.idx[k]])}
      }else{
        snp.val[,2]=.compute_snpQ2onepair(x@refallele.readcount[,pop.idx[2]],x@refallele.readcount[,pop.idx[3]],
                                          x@readcoverage[,pop.idx[2]],x@readcoverage[,pop.idx[3]])
        for(k in 2:3){#ab et ac
          snp.val[,k+1]=.compute_snpQ2onepair(x@refallele.readcount[,pop.idx[1]],x@refallele.readcount[,pop.idx[k]],
                                              x@readcoverage[,pop.idx[1]],x@readcoverage[,pop.idx[k]])}
      }
      snp.val=0.5*(rowSums(snp.val[,1:2]) - rowSums(snp.val[,3:4]))
    }    
    
    if(statname%in%c("F4","Dstat","Fh","Fd","FdM")){#all ahve the same numerator
      q2.req=cbind(pop.idx[c(1,2,1,2)],pop.idx[c(3,4,4,3)])#Q2ac Q2bd Q2ad Q2bc
      snp.val=matrix(0,x@nsnp,nrow(q2.req)) 
      for(i in 1:nrow(q2.req)){
        if(is.countdata(x)){
          snp.val[,i]=.compute_snpQ2onepair(x@refallele.count[,q2.req[i,1]],x@refallele.count[,q2.req[i,2]],
                                            x@total.count[,q2.req[i,1]],x@total.count[,q2.req[i,2]])
        }else{
          snp.val[,i]=.compute_snpQ2onepair(x@refallele.readcount[,q2.req[i,1]],x@refallele.readcount[,q2.req[i,2]],
                                            x@readcoverage[,q2.req[i,1]],x@readcoverage[,q2.req[i,2]])
        }
      }
      snp.val=0.5*(rowSums(snp.val[,1:2]) - rowSums(snp.val[,3:4]))
    }
    
    if(statname%in%c("DenFd","DenFdM")){#denominator is an F3 stats (x -1) which depends on allele freq
      tmp.F3conf=rbind(pop.idx[c(2,1,4)],pop.idx[c(3,1,4)],pop.idx[c(1,2,4)],pop.idx[c(3,2,4)]) 
      snp.val=rep(NA,x@nsnp)
      snp.class=rep(FALSE,x@nsnp)
      if(is.countdata(x)){
        tmp.freq=x@refallele.count[,pop.idx[1:4]]/x@total.count[,pop.idx[1:4]]
      }else{
        tmp.freq=x@refallele.readcount[,pop.idx[1:4]]/x@readcoverage[,pop.idx[1:4]]
      }
      #polarisation to avoid condition in determining auxiliary (i.e. allele ref = allele der (the less frequence in P4=Outgroup))
      tmp.revsnp=tmp.freq[,4]>0.5
      tmp.freq[tmp.revsnp,]=1-tmp.freq[tmp.revsnp,]
      rm(tmp.revsnp)
      if(statname=="DenFd"){
        tmp.nclass=2
        snp.class[rowSums(is.na(tmp.freq))>0]=NA
        snp.class[tmp.freq[,2]>tmp.freq[,3]]=1   #Fd=-F4(P1,P2;P3,P4)/F3(P2;P1,P4)
        snp.class[tmp.freq[,2]<=tmp.freq[,3]]=2  #Fd=-F4(P1,P2;P3,P4)/F3(P3;P1,P4)
      }else{
        tmp.nclass=4
        snp.class[rowSums(is.na(tmp.freq))>0]=NA
        snp.class[tmp.freq[,2]>=tmp.freq[,1] & tmp.freq[,2]>tmp.freq[,3]]=1  #FdM=Fd=-F4(P1,P2;P3,P4)/F3(P2;P1,P4)
        snp.class[tmp.freq[,2]>=tmp.freq[,1] & tmp.freq[,2]<=tmp.freq[,3]]=2 #FdM=Fd=-F4(P1,P2;P3,P4)/F3(P3;P1,P4)
        snp.class[tmp.freq[,2]<tmp.freq[,1] & tmp.freq[,1]>=tmp.freq[,3]]=3  #FdM=Fd=-F4(P1,P2;P3,P4)/F3(P1;P2,P4)
        snp.class[tmp.freq[,2]<tmp.freq[,1] & tmp.freq[,1]<tmp.freq[,3]]=4   #FdM=Fd=-F4(P1,P2;P3,P4)/F3(P3;P2,P4)
      }
      
      snp.val=matrix(NA,x@nsnp,4) #Q1a Q2bc Q2ab Q2ac
      for(tmp.class in 1:tmp.nclass){
        tmp.snp.idx=which(snp.class==tmp.class)
        if(length(tmp.snp.idx)>0){
          if(is.countdata(x)){
            #get info for Q1
            tmp.w=c(1)
            tmp.refcount=x@refallele.count[tmp.snp.idx,tmp.F3conf[tmp.class,1],drop=F]
            tmp.totcount=x@total.count[tmp.snp.idx,tmp.F3conf[tmp.class,1],drop=F]
            #computeQ2
            snp.val[tmp.snp.idx,2]=.compute_snpQ2onepair(x@refallele.count[tmp.snp.idx,tmp.F3conf[tmp.class,2]],
                                                         x@refallele.count[tmp.snp.idx,tmp.F3conf[tmp.class,3]],
                                                         x@total.count[tmp.snp.idx,tmp.F3conf[tmp.class,2]],
                                                         x@total.count[tmp.snp.idx,tmp.F3conf[tmp.class,3]])
            for(k in 2:3){#ab et ac
              snp.val[tmp.snp.idx,k+1]=.compute_snpQ2onepair(x@refallele.count[tmp.snp.idx,tmp.F3conf[tmp.class,1]],
                                                             x@refallele.count[tmp.snp.idx,tmp.F3conf[tmp.class,k]],
                                                             x@total.count[tmp.snp.idx,tmp.F3conf[tmp.class,1]],
                                                             x@total.count[tmp.snp.idx,tmp.F3conf[tmp.class,k]])}
          }else{
            #get info for Q1
            tmp.w=x@poolsizes[tmp.F3conf[tmp.class,1]]/(x@poolsizes[tmp.F3conf[tmp.class,1]]-1)
            tmp.refcount=x@refallele.readcount[tmp.snp.idx,tmp.F3conf[tmp.class,1],drop=F]
            tmp.totcount=x@readcoverage[tmp.snp.idx,tmp.F3conf[tmp.class,1],drop=F] 
            #computeQ2
            #computeQ2
            snp.val[tmp.snp.idx,2]=.compute_snpQ2onepair(x@refallele.readcount[tmp.snp.idx,tmp.F3conf[tmp.class,2]],
                                                         x@refallele.readcount[tmp.snp.idx,tmp.F3conf[tmp.class,3]],
                                                         x@readcoverage[tmp.snp.idx,tmp.F3conf[tmp.class,2]],
                                                         x@readcoverage[tmp.snp.idx,tmp.F3conf[tmp.class,3]])
            for(k in 2:3){#ab et ac
              snp.val[tmp.snp.idx,k+1]=.compute_snpQ2onepair(x@refallele.readcount[tmp.snp.idx,tmp.F3conf[tmp.class,1]],
                                                             x@refallele.readcount[tmp.snp.idx,tmp.F3conf[tmp.class,k]],
                                                             x@readcoverage[tmp.snp.idx,tmp.F3conf[tmp.class,1]],
                                                             x@readcoverage[tmp.snp.idx,tmp.F3conf[tmp.class,k]])}
          }    
          snp.val[tmp.snp.idx,1]=.compute_snpQ1onepop(tmp.refcount,tmp.totcount,tmp.w)
        }
      }
      snp.val= -0.5*(rowSums(snp.val[,1:2]) - rowSums(snp.val[,3:4]))
    }   
    
    
    return(snp.val)
  }
  
  ##########
  ##########
  
  if(verbose){cat("Computing SNP-specific values\n")} #needs to be done first (to account for NA in windows definition)
  if(verbose){cat("For the num.pop.idx combination\n")}
  snp.num.val=compute_snp_stat_(num.stat,num.pop.idx)
  keep=!is.na(snp.num.val)
  num.scaling=FALSE
  #denominator of the numerator! to properly scale stats when computing multilocus window stat
  if(num.stat%in%c("Fst","F3star","Dstat","Fh","Fd","FdM")){
    num.scaling=TRUE
    if(num.stat=="Fst"){snp.num.val.scale=compute_snp_stat_("div",num.pop.idx)}
    if(num.stat=="F3star"){snp.num.val.scale=compute_snp_stat_("het",num.pop.idx[1])}
    if(num.stat=="Dstat"){#/(1-Q2ab)*(1-Q2bc)
      snp.num.val.scale=compute_snp_stat_("div",num.pop.idx[1:2])*compute_snp_stat_("div",num.pop.idx[3:4])}
    if(num.stat=="Fh"){snp.num.val.scale=-1*compute_snp_stat_("F3",num.pop.idx[c(3,1,4)])}
    if(num.stat=="Fd"){snp.num.val.scale=compute_snp_stat_("DenFd",num.pop.idx)} # *(-1) done in compute_sno_stat_
    if(num.stat=="FdM"){snp.num.val.scale=compute_snp_stat_("DenFdM",num.pop.idx)} # *(-1) done in compute_sno_stat_
    keep=keep & !is.na(snp.num.val.scale)
  }
  
  if(denom){
    if(verbose){cat("For the den.pop.idx combination\n")}
    snp.den.val=compute_snp_stat_(den.stat,den.pop.idx)
    keep=keep & !is.na(snp.den.val)
    den.scaling=FALSE
    if(den.stat%in%c("Fst","F3star","Dstat","Fh","Fd","FdM")){
      den.scaling=TRUE
      if(den.stat=="Fst"){snp.den.val.scale=compute_snp_stat_("div",den.pop.idx)}
      if(den.stat=="F3star"){snp.den.val.scale=compute_snp_stat_("het",den.pop.idx[1])}
      if(den.stat=="Dstat"){#/(1-Q2ab)*(1-Q2bc)
        snp.den.val.scale=compute_snp_stat_("div",den.pop.idx[1:2])*compute_snp_stat_("div",den.pop.idx[3:4])}
      if(den.stat=="Fh"){snp.den.val.scale=-1*compute_snp_stat_("F3",den.pop.idx[c(3,1,4)])}
      if(den.stat=="Fd"){snp.den.val.scale=compute_snp_stat_("DenFd",den.pop.idx)} # *(-1) done in compute_snp_stat_
      if(den.stat=="FdM"){snp.den.val.scale=compute_snp_stat_("DenFdM",den.pop.idx)} # *(-1) done in compute_snp_stat_
      keep=keep & !is.na(snp.den.val.scale)
      snp.den.val.scale=snp.den.val.scale[keep]
    }
    snp.den.val=snp.den.val[keep]
  }
  snp.num.val=snp.num.val[keep]
  if(num.scaling){snp.num.val.scale=snp.num.val.scale[keep]}
  
  ########  
  snp.keep.info=x@snp.info[keep,]
  if(verbose){cat("Defining Windows\n")}
  #Window definition
  win.det=matrix(NA,0,9)
  colnames(win.det)=c("Chr","Start","End","MidPos","CumMidPos","id_snp_start","id_snp_end","nsnp","Value")
  win.det=as.data.frame(win.det)
  chr=unique(as.character(snp.keep.info$Chromosome))
  nchr=length(chr)
  cum.midpos=0
  for(c in 1:nchr){
    tmp.continue=TRUE
    tmp.sel=which(snp.keep.info$Chromosome==chr[c])
    tmp.pos=floor(snp.keep.info$Position[tmp.sel]) #in case real (e.g., when simulated data)
    tmp.nsnps=length(tmp.sel)
    if(winbp){
      if(bp.start.first.snp){
        pos.init=tmp.pos[1]
      }else{
        pos.init=(tmp.pos[1]%/%sliding.window.size)*sliding.window.size + 1
        if(tmp.pos[1]%%sliding.window.size==0){pos.init=pos.init-sliding.window.size}#the first SNP is the only of the first window (that may be discarded further!)
      }
      if((tmp.pos[tmp.nsnps]-pos.init)>sliding.window.size){
        tmp.start.pos=seq(pos.init,tmp.pos[tmp.nsnps],step)
        tmp.end.pos = tmp.start.pos + sliding.window.size
        tmp.snp.start.idx=sapply(tmp.start.pos,f<-function(x){min(which(tmp.pos>=x))})
        tmp.snp.end.idx=sapply(tmp.end.pos,f<-function(x){max(which(tmp.pos<=x))}) 
      }else{
        tmp.continue=FALSE
      }
    }else{
      if(tmp.nsnps>sliding.window.size){
        tmp.snp.start.idx=seq(1,tmp.nsnps,step)
        tmp.snp.end.idx=tmp.snp.start.idx+sliding.window.size-1
        tmp.snp.end.idx=pmin(tmp.snp.end.idx,tmp.nsnps)
        tmp.start.pos=tmp.pos[tmp.snp.start.idx]
        tmp.end.pos=tmp.pos[tmp.snp.end.idx]
      }else{
        tmp.continue=FALSE
      }
    }
    if(tmp.continue){
      tmp.mid.pos=(tmp.start.pos+tmp.end.pos)/2
      tmp.cum.mid.pos=tmp.mid.pos+cum.midpos
      cum.midpos=cum.midpos+max(tmp.mid.pos)
      tmp.win.nsnp=tmp.snp.end.idx-tmp.snp.start.idx+1
      tmp.nwin=length(tmp.start.pos)
      tmp.windet=data.frame(Chr=rep(chr[c],tmp.nwin),Start=tmp.start.pos,End=tmp.end.pos,
                            MidPos=tmp.mid.pos,CumMidPos=tmp.cum.mid.pos,
                            id_snp_start=tmp.sel[tmp.snp.start.idx],
                            id_snp_end=tmp.sel[tmp.snp.end.idx],
                            nsnp=tmp.win.nsnp,Value=rep(NA,tmp.nwin))
      if(winbp){
        tmp.win.size=tmp.windet$End-tmp.windet$Start
        tmp.win.keep=(tmp.win.size>0.99*sliding.window.size) & (tmp.win.size<1.01*sliding.window.size) & 
          tmp.windet$nsnp>1
      }else{
        tmp.win.keep=tmp.windet$nsnp==sliding.window.size
      }
      win.det=rbind(win.det,tmp.windet[tmp.win.keep,])               
    }
  }
  
  nwins=nrow(win.det)                         
  if(nwins<1){
    stop("No SNP-window found with the given parameters. Try decreasing sliding.window.size\n")
  }
  
  if(verbose){
    cat(nwins,"(overlapping) windows identified\n")
    dum=(win.det$End-win.det$Start)*1e-3
    cat("  Average (min-max) window sizes (in kb):",round(mean(dum),1),"(",round(min(dum)),"-",round(max(dum)),")\n")
    cat("  Average (min-max) nb. of SNPs per window:",round(mean(win.det$nsnp),1),"(",min(win.det$nsnp),"-",max(win.det$nsnp),")\n")    
  }
  ####### 
  
  if(verbose){cat("Computing window statistics\n")}
  snp.win.idx=cbind(win.det$id_snp_start,win.det$id_snp_end)-1 #0-indexed (for cpp)
  win.det$Value=.block_sum2(snp.num.val,snp.win.idx)
  if(num.scaling){win.det$Value=win.det$Value/.block_sum2(snp.num.val.scale,snp.win.idx)}
  if(denom){
    win.det$Value=win.det$Value/.block_sum2(snp.den.val,snp.win.idx)
    if(den.scaling){win.det$Value=win.det$Value * .block_sum2(snp.den.val.scale,snp.win.idx)}
  }else{
    if(!num.scaling){win.det$Value=win.det$Value/win.det$nsnp}
  }
  return(win.det[,c(1:5,8:9)])
}