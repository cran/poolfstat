#' Find sets of populations that may used as scaffold tree
#' @param fstats Object of class fstats containing estimates of fstats (see the function compute.fstats)
#' @param f3.zcore.threshold The significance threshold for Z-score of formal test of admixture based on the F3-statistics (default=-2)
#' @param f4.zscore.absolute.threshold The significance threshold for |Z-score| of formal test of treeness based on the F4-statistics  (default=2)
#' @param excluded.pops Vector of pop names to be exclude from the exploration
#' @param nthreads Number of available threads for parallelization of some part of the parsing (default=1, i.e., no parallelization)
#' @param verbose If TRUE extra information is printed on the terminal
#' @details The procedure first discards all the populations P that shows a significant signal of admixture with a Z-score for F3 statistics of the form F3(P;Q,R) < f3.zscore.thresholds. It then identifies all the sets of populations that pass the F4-based treeness with themselves. More precisely, for a given set E containing n populations, the procedure ensure that all the n(n-1)(n-2)(n-3)/8 possible F4 quadruplets have a |Z-score|<f4.zscore.absolute.threshold. The function aims at maximizing the size of the sets. 
#' @return A list with the following elements:
#' \enumerate{
#' \item "n.sets": The number of sets of (scaffold) unadmixed populations identified
#' \item "set.size": The number of populations included in each set
#' \item "pop.sets": A character matrix of n.sets rows and set.size columns giving for each set identified the names of the included populations.
#' \item "Z_f4.range": A matrix of n.sets rows and 2 columns reported for each set the range of variation (min and max value) of the absolute F4 Z-scores for the quadruplets passing the treeness test. More precisely, for a given set consisting of n=set.size populations, a total of n(n-1)(n-2)(n-3)/8 quadruplets can be formed. Yet, any set of four populations A, B, C and D is represented by three quadruplets A,B;C,D (or one of its seven other equivalent combinations formed by permuting each pairs); A,C;B,D (or one of its seven other equivalent combinations) and A,D;B,C (or one of its seven other combinations). Among these three, only a single quadruplet is expected to pass the treeness test (i.e., if the correct unrooted tree topology is (A,C;B,D), then the absoulte value of the Z-scores associated to F4(A,B;C,D) and F4(A,D;B,C) or their equivalent will be high. 
#' \item "passing.quadruplets": A matrix of n.sets rows and set.size columns reporting for each sets the n(n-1)(n-2)(n-3)/24 quadruplets that pass the treeness test (see Z_f4.range detail). 
#' }
#' @seealso see \code{\link{compute.fstats}}.
#' @examples
#' make.example.files(writing.dir=tempdir())
#' pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#' res.fstats=compute.fstats(pooldata,nsnp.per.bjack.block = 50) 
#' #NOTE: toy example (in practice nsnp.per.bjack.block should be higher)
#' popsets=find.tree.popset(res.fstats,f3.zcore.threshold=-3)  
#' @export
find.tree.popset<-function(fstats,f3.zcore.threshold=-1.65,f4.zscore.absolute.threshold=1.96,excluded.pops=NULL,nthreads=1,verbose=TRUE){
  if(nthreads>1){
    tmp.ncores=detectCores()
    if(nthreads>tmp.ncores){nthreads=tmp.ncores}
    options(cores=nthreads)
    registerDoParallel()  ;  getDoParWorkers()  
    parallel=TRUE
  }else{parallel=FALSE}
  if(!(is.fstats(fstats))){stop("Input data should either be an object of class fstats (see the function compute.fstats)\n")}
  if(!fstats@blockjacknife){stop("fstats object must contain estimates of Z-score for F3 and F4-based criteria to be evaluated (see compute.fstats taht should be run with nsnp.per.bjack.block>0 to allow block-jackknife estimates of s.e. of the estimates)\n")}
  if(nrow(fstats@f3.values)==0){stop("fstats object must contain estimates of F3 (see compute.fstats that should be run with computeF3=TRUE)\n")}
  if(nrow(fstats@f4.values)==0){stop("fstats object must contain estimates of F4 (see compute.fstats that should be run with computeF4=TRUE)\n")}

  f3.pop.names=fstats@comparisons[["F3"]]
  f4.pop.names=fstats@comparisons[["F4"]] 
  
  #Find quartets passing F4 treenes test
  comp.kept=abs(fstats@f4.values[,"Z-score"])<f4.zscore.absolute.threshold
  if(sum(comp.kept)==0){
    stop("No quadruplets could be found (none passing the F4 treeness test)\n")
  }
  if(sum(comp.kept)==1){#otherwise vector and problems in the remaining
    pop.quadruplets=t(matrix(f4.pop.names[comp.kept,]))
  }else{
   pop.quadruplets=f4.pop.names[comp.kept,]
  }
  #exclude pops in excluded.pops vector
  if(!is.null(excluded.pops)){
    comp.kept=apply(pop.quadruplets,1,f<-function(x){sum(x %in% excluded.pops)==0})
    if(sum(comp.kept)==0){
      stop("No quadruplets left after excluding those involving populations in the vector excluded.pops\n")
    }
    if(sum(comp.kept)==1){#otherwise vector and problems in the remaining
      pop.quadruplets=t(matrix(pop.quadruplets[comp.kept,]))
    }else{
      pop.quadruplets=pop.quadruplets[comp.kept,]
    }
  }
  
  #Remove those containing a pop that pass the F3 tests (if any)
  #Find positive F3 tests (population will be automatically excluded)
  comp.elim=fstats@f3.values[,"Z-score"]<f3.zcore.threshold
  if(sum(comp.elim)>0){
    pop.elim=unique(f3.pop.names[comp.elim,1])
    comp.kept=apply(pop.quadruplets,1,f<-function(x){sum(x %in% pop.elim)==0})
    if(sum(comp.kept)==0){
      stop("No quadruplets remaining after F3 filtering\n")
    }else{
      if(sum(comp.kept)==1){#otherwise vector and problems in the remaining
        pop.quadruplets=t(matrix(pop.quadruplets[comp.kept,]))
      }else{
        pop.quadruplets=pop.quadruplets[comp.kept,]
      }
    }
  }

  n.quads=nrow(pop.quadruplets)
  if(n.quads==1){#then a asingle qusruplet remaining at this stage!
    subset.pops=pop.quadruplets
  }else{
   all.pops=unique(as.character(pop.quadruplets))
   #incidence matrix (helpful for testing F4 tests)
   mat.f4.ok=t(apply(pop.quadruplets,1,f<-function(x){all.pops %in% x}))
   colnames(mat.f4.ok)=all.pops
   subset.new=pop.quadruplets
   while(nrow(subset.new)>0 & ncol(subset.new)<length(all.pops)){
    if(verbose){cat("Number of sets:",nrow(subset.new),"of Npops=",ncol(subset.new),"each\n")}
    all.pops=unique(as.character(subset.new))
    subset.cur=subset.new
    subset.length=ncol(subset.cur)
    if(parallel){
      subset.new=foreach(i=1:nrow(subset.cur),.combine=c) %dopar% {      
        tmp.subset=subset.cur[i,]
        tmp.pop.test=all.pops[!(all.pops %in% tmp.subset)]
        tmp.subset.triplets=combn(tmp.subset,3)
        tmp.ntriplets=ncol(tmp.subset.triplets)
        tmp.subset.new=c()
        for(j in tmp.pop.test){
          tmp.cnt=0
          for(k in 1:ncol(tmp.subset.triplets)){
            if(sum(rowSums(mat.f4.ok[,c(tmp.subset.triplets[,k],j)])==4)>0){tmp.cnt=tmp.cnt+1}
          }
          if(tmp.cnt==tmp.ntriplets){tmp.subset.new=c(tmp.subset.new,paste(sort(c(tmp.subset,j)),collapse=":"))}
        } 
        tmp.subset.new
      }
    }else{
      subset.new=c()    
      for(i in 1:nrow(subset.cur)){
        tmp.subset=subset.cur[i,]
        tmp.pop.test=all.pops[!(all.pops %in% tmp.subset)]
        tmp.subset.triplets=combn(tmp.subset,3)
        tmp.ntriplets=ncol(tmp.subset.triplets)
        for(j in tmp.pop.test){
          tmp.cnt=0
          for(k in 1:ncol(tmp.subset.triplets)){
            if(sum(rowSums(mat.f4.ok[,c(tmp.subset.triplets[,k],j)])==4)>0){tmp.cnt=tmp.cnt+1}
          }
          if(tmp.cnt==tmp.ntriplets){subset.new=c(subset.new,paste(sort(c(tmp.subset,j)),collapse=":"))}
        }
      }
    }
    if(length(subset.new)>0){
      subset.new=unique(subset.new)
      subset.new=matrix(unlist(strsplit(subset.new, split=":")),ncol=subset.length+1,byrow=T)
    }else{
      subset.new=matrix("",nrow=0,ncol=subset.length+1)
    }
  }
  if(nrow(subset.new)==0){subset.pops=subset.cur}else{subset.pops=subset.new}
  }
  ####
  #return stats for the subset.pops (min and max F4 Z-scores)
  n.subsets=nrow(subset.pops)
  summary.f4.zscores=matrix(0,n.subsets,2)
  colnames(summary.f4.zscores)=c("Min. |Zscore|","Max. |Zscore|")
  for(i in 1:n.subsets){
    tmp.comp = apply(fstats@comparisons[["F4"]] ,1,f<-function(x){sum(x %in% subset.pops[i,])==4})
    tmp.f4zscores=fstats@f4.values$`Z-score`[tmp.comp]
    #for every quadruplets only 1 among the 3 must pass the treeness test: if (A,B;C,D) is the topology OK then (A,C;B,D) (and the seven other combinations) and (A,D;B,C) (and the 7 other combinations) must have a very high Z-scores. Instead of exploring we simply take the n.zscores/3 lowest zscores to compute the summary
    tmp.conf.id=order(abs(tmp.f4zscores),decreasing = F)[1:(length(tmp.f4zscores)/3)]
    summary.f4.zscores[i,]=range(abs(tmp.f4zscores[tmp.conf.id]))
    if(i==1){
      conf.ok=rownames(fstats@comparisons[["F4"]])[tmp.comp][tmp.conf.id]
    }else{
      conf.ok=rbind(conf.ok,rownames(fstats@comparisons[["F4"]])[tmp.comp][tmp.conf.id])
    }
  }
  
  if(n.subsets==1){conf.ok=t(as.matrix(conf.ok))}
  rownames(conf.ok)=rownames(subset.pops)=rownames(summary.f4.zscores)=paste0("PopSet",1:n.subsets)
  
  return(list(n.sets=n.subsets,set.size=ncol(subset.pops),pop.sets=subset.pops,Z_f4.range=summary.f4.zscores,passing.quaduplets=conf.ok))
}

