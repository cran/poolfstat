#' Compare fitted f2, f3 and f4 f-statistics of an admixture graph with estimated ones
#' @param fstats Object of class fstats containing estimates of fstats (as obtained with compute.fstats)
#' @param fitted.graph Object of class fitted graph (as obtained with fit.graph function).
#' @param n.worst.stats The number of worst statistics to be displayed in the terminal 
#' @details Compare fitted and estimated f-statistics may allow identifying problematic edges on the graph.   
#' @return A matrix with 3 columns for each test (row names of the matrix corresponding to the test):
#' \enumerate{
#'  \item The estimated f-statistics (mean across block-Jackknife samples)
#'  \item The fitted f-statistics (obtained from the fitted grah parameters
#'  \item A Z-score measuring the deviation of the fitted values from the estimated values in units of standard errors (i.e., Z=(fitted.value-target.value)/se(target.value))
#' }
#' @seealso See \code{\link{compute.fstats}} and \code{\link{fit.graph}}
#' @export
compare.fitted.fstats<-function(fstats,fitted.graph,n.worst.stats=5){
  if(!is.fitted.graph(fitted.graph)){stop("The input fitted.graph is not a valid fitted.graph object (see the function fit.graph)\n")}
  if(!is.fstats(fstats)){stop("The input fstats is not a valid fstats object (see the function compute.fstats)\n")}  
  
  if(nrow(fstats@f3.values)==0){stop("No estimates of F3 statistics available in the fstats object.\nFunctions compute.fstats must be run with computeF3=TRUE\n")}
  if(nrow(fstats@f4.values)==0){stop("No estimates of F4 statistics available in the fstats object.\nFunctions compute.fstats must be run with computeF4=TRUE\n")}  

  npop=nrow(fitted.graph@fitted.f2.mat)
  pop.names=rownames(fitted.graph@fitted.f2.mat)
  pop.idx=1:npop
  names(pop.idx)=pop.names
  n.f2.exp=npop*(npop-1)/2
  n.f3.exp=npop*(npop-1)*(npop-2)/2
  n.f4.exp=((npop*(npop-1)/2) * ((npop-2)*(npop-3)/2)) / 2 
  
  # f2.mat=diag(fitted.graph@fitted.omega)%*%t(rep(1,nrow(fitted.graph@fitted.omega)))
  # f2.mat=f2.mat+t(f2.mat)-2*fitted.graph@fitted.omega
  f2.mat=fitted.graph@fitted.f2.mat
  
  f2.pops=fstats@comparisons[["F2"]]
  tst.sel=apply(f2.pops,1,f<-function(x){sum(x %in% pop.names)==2})
  if(sum(tst.sel)!=n.f2.exp){
    cat(n.f2.exp,"F2 expected and only",sum(tst.sel),"given in the fstats.estimated object\n")
    warning("Some F2 configurations are missing\n")
    if(sum(tst.sel)==0){stop("ERROR: No F2 in common between fitted and estimated objects\n")}
  }
  f2.pops=f2.pops[tst.sel,]
  f2.idx=(pop.idx[f2.pops[,1]]-1)*npop+pop.idx[f2.pops[,2]]
  f2.fitted=f2.mat[f2.idx]
  f2.fitted=cbind(fstats@f2.values[tst.sel,2],f2.fitted,(f2.fitted-fstats@f2.values[tst.sel,2])/fstats@f2.values[tst.sel,3])
  rownames(f2.fitted)=rownames(fstats@f2.values[tst.sel,])
    
  f3.pops=fstats@comparisons[["F3"]]
  tst.sel=apply(f3.pops,1,f<-function(x){sum(x %in% pop.names)==3})
  if(sum(tst.sel)!=n.f3.exp){
    cat(n.f3.exp,"F3 expected and only",sum(tst.sel),"given in the fstats.estimated object\n")
    warning("Some F3 configurations are missing\n")
    if(sum(tst.sel)==0){stop("ERROR: No F3 in common between fitted and estimated objects\n")}
  }
  f3.pops=f3.pops[tst.sel,]
  f2.for.f3.idx.1=(pop.idx[f3.pops[,1]]-1)*npop+pop.idx[f3.pops[,2]]
  f2.for.f3.idx.2=(pop.idx[f3.pops[,1]]-1)*npop+pop.idx[f3.pops[,3]]
  f2.for.f3.idx.3=(pop.idx[f3.pops[,2]]-1)*npop+pop.idx[f3.pops[,3]]
  f3.fitted=(f2.mat[f2.for.f3.idx.1]+f2.mat[f2.for.f3.idx.2]-f2.mat[f2.for.f3.idx.3])/2
  f3.fitted=cbind(fstats@f3.values[tst.sel,2],f3.fitted,(f3.fitted-fstats@f3.values[tst.sel,2])/fstats@f3.values[tst.sel,3])
  rownames(f3.fitted)=rownames(fstats@f3.values[tst.sel,])
  
  f4.pops=fstats@comparisons[["F4"]]
  tst.sel=apply(f4.pops,1,f<-function(x){sum(x %in% pop.names)==4})
  if(sum(tst.sel)!=n.f4.exp){
    cat(n.f4.exp,"F4 expected and only",sum(tst.sel),"given in the fstats.estimated object\n")
    warning("Some F4 configurations are missing\n")
    if(sum(tst.sel)==0){stop("ERROR: No F4 in common between fitted and estimated objects\n")}
  }
  f4.pops=f4.pops[tst.sel,]
  f2.for.f4.idx.1=(pop.idx[f4.pops[,1]]-1)*npop+pop.idx[f4.pops[,3]]
  f2.for.f4.idx.2=(pop.idx[f4.pops[,1]]-1)*npop+pop.idx[f4.pops[,4]]
  f2.for.f4.idx.3=(pop.idx[f4.pops[,2]]-1)*npop+pop.idx[f4.pops[,3]]
  f2.for.f4.idx.4=(pop.idx[f4.pops[,2]]-1)*npop+pop.idx[f4.pops[,4]]
  f4.fitted=(f2.mat[f2.for.f4.idx.2]+f2.mat[f2.for.f4.idx.3]-f2.mat[f2.for.f4.idx.1]-f2.mat[f2.for.f4.idx.4])/2
  f4.fitted=cbind(fstats@f4.values[tst.sel,2],f4.fitted,(f4.fitted-fstats@f4.values[tst.sel,2])/fstats@f4.values[tst.sel,3])
  rownames(f4.fitted)=rownames(fstats@f4.values[tst.sel,])
    
  comp.res=rbind(f2.fitted,f3.fitted,f4.fitted)
  colnames(comp.res)=c("Estimated","Fitted","Z--score")
  comp.res=as.data.frame(comp.res)
  if(n.worst.stats>0){
   n.worst.stats=min(max(1,n.worst.stats),nrow(comp.res))
   cat(n.worst.stats,"Worst fit for:\n")
   show(comp.res[order(abs(comp.res[,3]),decreasing = T)[1:n.worst.stats],])
  }
  return(comp.res)
}
