#' Plot F2, F3, F3star, F4, D or pairwise Fst values with their Confidence Intervals
#' @param x An object of class fstats (to plot F2, F3 or F4 statistics) or pairwisefst (to plot pairwise fst)
#' @param stat.name For fstats object, the name of the stat (either F2, F3, F3star, F4 or Dstat)
#' @param ci.perc Percentage of the Confidence Interval in number of standard errors (default=95\%)
#' @param value.range Range of test values (x-axis) to be plotted (default=NA,NA: i.e., all test values are plotted)
#' @param pop.sel Only plot test values involving these populations (default=NA: i.e., all test values are plotted)
#' @param pop.f3.target For F3-statistics, only plot F3 involving pop.f3.target as a target
#' @param highlight.signif If TRUE highlight significant tests in red (see details)
#' @param main Main title of the plot (default=stat.name)
#' @param ... Some other graphical arguments to be passed
#' @details Data will only be plotted if jackknife estimates of the estimator s.e. have been performed i.e. if the functions compute.fstats or compute.pairwiseFST were run with nsnp.per.block>0
#' @return A plot of the Fstats of interest. Significant F3 statistics (i.e., showing formal evidence for admixture of the target population) are highlighted in red. Significant F4 statistics (i.e., showing formal evidence against treeness of the pop. quadruplet) are highlighted in red. 
#' @seealso To generate x object, see \code{\link{compute.pairwiseFST}} (for pairwisefst object) or \code{\link{compute.fstats}}  (for fstats object)
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),
#'                    poolsizes=rep(50,15),poolnames=paste0("P",1:15))
#'  res.fstats=compute.fstats(pooldata,nsnp.per.bjack.block=25)
#'  plot_fstats(res.fstats,stat.name="F3",cex=0.5)
#'  plot_fstats(res.fstats,stat.name="F3",value.range=c(NA,0.001),
#'              pop.f3.target=c("P7","P5"),cex.axis=0.7)
#'  plot_fstats(res.fstats,stat.name="F4",cex=0.5) 
#'  #allow to reduce the size of the test name (y-axis)
#'  plot_fstats(res.fstats,stat.name="F4",cex=0.5,
#'              pop.sel=c("P1","P2","P3","P4","P5")) 
#'  plot_fstats(res.fstats,stat.name="F4",cex=0.5,
#'              pop.sel=c("P1","P2","P3","P4","P5"),highlight.signif=FALSE)  
#' @export
plot_fstats <- function(x, stat.name="F2",ci.perc = 95,value.range=c(NA,NA),pop.sel=NA,pop.f3.target=NA,highlight.signif=TRUE,main=stat.name, ... ) {
  if(!(is.pairwisefst(x)) & !(is.fstats(x))){
    stop("Input data should either be an object of class pairwisefst (see the function compute.pairwiseFST) or of class fstats (see the function compute.fstats)\n")} 
  if(!x@blockjacknife){
    stop("No block-jackknife estimate of estimators s.e.. Functions compute.fstats or compute.pairwiseFST should be run with nsnp.per.block>0\n")
  }
  
  if(is.pairwisefst(x)){
    tmp.x=x@values[,2:3]
    dum.pops=matrix(unlist(strsplit(rownames(tmp.x),split="[;]")),nrow=nrow(tmp.x),byrow=T)
    stat.name="Fst"
  }else{
      if(!(stat.name %in% c("F2","F3","F4","Fst","F3star","Dstat"))){
        cat("stat.name should either be equal to F2, F3, F3star, F4, F4star or Fst for a fstats object to be plotted\n")
        stop("")}
    if(stat.name=="F2"){tmp.x=x@f2.values[,2:3]}
    if(stat.name=="Fst"){tmp.x=x@fst.values[,2:3]}    
    if(stat.name=="F3"){tmp.x=x@f3.values[,2:3]}
    if(stat.name=="F3star"){tmp.x=x@f3star.values[,2:3]}    
    if(stat.name=="F4"){tmp.x=x@f4.values[,2:3]}
    if(stat.name=="Dstat"){tmp.x=x@Dstat.values[,2:3]}
    dum.pops=x@comparisons[[stat.name]]
    }
  
  if(!(ci.perc>50 & ci.perc<100)){
    warning("CI percentage (argument ci.perc) must be between 50 and 100 (set to default value of 95)\n")
    ci.size=95
  }
  ci.size=qnorm((1+0.01*ci.perc)/2)
  tmp.stderr=ci.size*tmp.x[,2]
  tmp.data=cbind(tmp.x[,1],tmp.x[,1]-tmp.stderr,tmp.x[,1]+tmp.stderr)
  rownames(tmp.data)=rownames(tmp.x)
  tmp.ord=order(tmp.data[,1],decreasing=T)
  tmp.data=tmp.data[tmp.ord,] ;  dum.pops=dum.pops[tmp.ord,]
  value.range[1]=max(value.range[1],min(tmp.data[,2]),na.rm=T)
  value.range[2]=min(value.range[2],max(tmp.data[,3]),na.rm=T)  
  tmp.sel=tmp.data[,2]>=value.range[1]&tmp.data[,3]<=value.range[2]
  tmp.data=tmp.data[tmp.sel,]
  dum.pops=dum.pops[tmp.sel,] #required if additional selection

  #####Pop selection if any
  pop.lst=unique(as.character(dum.pops))
  n.pops.per.test=ncol(dum.pops)
  if(sum(is.na(pop.sel))==0){
    if(sum(!(pop.sel%in%pop.lst))>0){stop("ERROR: some pop names in pop.sel are not in the test names\n")}
    npopsel=length(pop.sel)
    tst.sel=apply(dum.pops,1,f<-function(z){sum(z%in%pop.sel)>=min(npopsel,n.pops.per.test)})
    if(sum(tst.sel)<1){stop("ERROR: no test found with the required population configuration (defined in pop.sel)\n")}
    tmp.data=tmp.data[tst.sel,]
    dum.pops=dum.pops[tst.sel,] #required if additional selection
  }
  ####
  if(sum(is.na(pop.f3.target))==0){
    if(!(stat.name%in%c("F3","F3star"))){
      warning("Option pop.f3.target disregarded since only compatible with F3 or F3star statistics")
    }else{
      if(sum(!(pop.f3.target%in%pop.lst))>0){stop("ERROR: some pop names in pop.f3.target are not in the test names\n")}
      tst.sel=dum.pops[,1]%in%pop.f3.target
      if(sum(tst.sel)<1){stop("ERROR: no test found with the required population configuration for F3 or F3star (defined in pop.f3.target)\n")}
      tmp.data=tmp.data[tst.sel,]    
      dum.pops=dum.pops[tst.sel,] #required if additional selection
    }
  }
  ####
  n.test=nrow(tmp.data)
  if(n.test<2){stop("ERROR: not enough stats to plot (check value.range)\n")}
  op=par()
  par(mar=c(3,10,1,1)+0.1)
  plot(tmp.data[,1],1:n.test,xlim=c(min(tmp.data[,2]),max(tmp.data[,3])),main=main,yaxt="n",ylab="",xlab="",pch="",bty="n",...)
  tmp.col=rep("black",n.test)
  if(highlight.signif){
   if(stat.name%in%c("F3","F3star")){tmp.col[tmp.data[,3]<0]="red"}
   if(stat.name%in%c("F4","Dstat")){tmp.col[tmp.data[,2]*tmp.data[,3]>0]="red"}  
  }
  text(labels=rownames(tmp.data), col=tmp.col, x=rep(min(tmp.data[,2]),n.test), y=1:n.test, pos=2,srt = 0, xpd = NA,...)
  sapply(1:n.test,f<-function(z){
    abline(h=z,lty=3,col="grey")
    arrows(tmp.data[z,2],z,tmp.data[z,3],z,angle = 90,code=3,length=0.01,col=tmp.col[z])
  })
  sapply(1:n.test,f<-function(z){abline(h=z,lty=3,col="grey")})  
  par(mar=op$mar)
  abline(v=0,lty=2,col="red")
}
