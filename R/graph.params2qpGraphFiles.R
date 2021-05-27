#' Generate files for the qpGraph software from a graph.params object
#' @param graph.params An object of class graph.params containing graph information with Fstats information (see the function generate.graph.params)
#' @param outfileprefix The prefix of the qpGraph files
#' @param n.printed.dec Number of decimal to be printed (if not enough may lead to fatalx error in qpGraph)
#' @param verbose If TRUE extra information is printed on the terminal
#' @details This function generates the three files required by qpGraph: i) a file named {outfileprefix}.graph containing the graph in appropriate format; ii) a file named {outfileprefix}.fstats file containing the fstats estimates of fstats (and their covariance); iii) a file named {outfileprefix}.parqpGraph containing essential parameter information to run qpGraph (this may be edited by hand if other options are needed). The qpGraph software may then be run using the following options -p {outfileprefix}.parqpGraph -g {outfileprefix}.graph -o out.ggg -d out.dot. 
#' @return The three files described in the details section
#' @seealso To generate graph.params object, see \code{\link{generate.graph.params}}
#' @export
graph.params2qpGraphFiles<-function(graph.params,outfileprefix="out",n.printed.dec=4,verbose=TRUE){
 ##Some checks
  if(!(is.graph.params(graph.params))){
    stop("The input graph.params is not a valid graph.params object (see generate.graph.params)\n")
  }
  if(length(graph.params@f2.target)==0){
    stop("The input graph.params does not contain fstats estimates (the function generate.graph.params to create it may have been run without fstats object)\n")    
  }
  if(!(n.printed.dec %in% 1:8)){stop("n.printed.dec must be an integer >=1 and <=8\n")}
  #to control "fatalx: (vlog): negative or zero value 0)" that may appear with qpGraph if cov two small
  f.prec=paste0("%10.",n.printed.dec,"f")
  covfact=1e6;f3fact=1e3
  if(nchar(outfileprefix)==0){outprefix="out"}
  #######################
  ####out.fstats file
  ######################
  if(is.null(graph.params@f2.target)){stop("Invalid graph.params object: see generate.graph.params function\n")}
  outfile=paste0(outfileprefix,".fstats")
  cat(file=outfile,paste0("##fbasis.  basepop:",graph.params@popref,"::  f3*",format(f3fact,scientific = F)," covar*",format(covfact,scientific = F),"\n"))
  for(i in 1:length(graph.params@f2.target)){
      cat(sprintf("%15s",graph.params@f2.target.pops[i,2]),
          sprintf("%15s",graph.params@f2.target.pops[i,2]),
          sprintf(f.prec,graph.params@f2.target[i]*f3fact),"\n",file=outfile,append=T)
  }
  for(i in 1:length(graph.params@f3.target)){
    cat(sprintf("%15s",graph.params@f3.target.pops[i,2]),
        sprintf("%15s",graph.params@f3.target.pops[i,3]),
        sprintf(f.prec,graph.params@f3.target[i]*f3fact),"\n",file=outfile,append=T)
  }
  tmp.n=length(graph.params@f2.target)+length(graph.params@f3.target)
  tmp.nomcov=rbind(cbind(graph.params@f2.target.pops[,2],graph.params@f2.target.pops[,2]),
                   graph.params@f3.target.pops[,2:3])
    for(i in 1:tmp.n){
      for(j in i:tmp.n){
        cat(sprintf("%15s",tmp.nomcov[i,1]),sprintf("%15s",tmp.nomcov[i,2]),
            sprintf("%15s",tmp.nomcov[j,1]),sprintf("%15s",tmp.nomcov[j,2]),
            sprintf(f.prec,graph.params@f.Qmat[i,j]*covfact),"\n",file=outfile,append=T)        
      }
    }
    if(verbose){cat("Fstats input file for qpGraph written in",outfile,"\n") } 
    #######################
    ####prepare graph file in qpGraph format
    ###################### 
    outgraphfile=paste0(outfileprefix,".graph")
    cat(file=outgraphfile,paste0("root\t",graph.params@root.name,"\n"))
    #for whatever reason: leaves should be ordered this way (otherwise error in qpGraph: "fatalx: (vlog): negative or zero value 0)
    cat(file=outgraphfile,paste0("label\t",graph.params@popref,"\t",graph.params@popref,"\n"),append=TRUE)
    for(i in 1:nrow(graph.params@f2.target.pops)){
      dum.pops=graph.params@f2.target.pops[i,2]
      cat(file=outgraphfile,paste0("label\t",dum.pops,"\t",dum.pops,"\n"),append=TRUE)
    }
    cat(file=outgraphfile,"\n",append=TRUE)
    tmp.graph=graph.params@graph
    adm.graph.rows=which(nchar(tmp.graph[,3])>0)
    if(length(adm.graph.rows)>0){
      adm.pops=unique(tmp.graph[adm.graph.rows,1])
      for(i in adm.pops){
       tmp.par=tmp.graph[tmp.graph[,1]==i,2]
       cat(file=outgraphfile,paste0("admix\t",i,"\t",tmp.par[1],"\t",tmp.par[2],"\n"),append=TRUE)
      }
      tmp.graph=tmp.graph[-1*adm.graph.rows,]  
    }    
    for(i in 1:nrow(tmp.graph)){
      cat(file=outgraphfile,paste0("edge\t",paste(tmp.graph[i,2:1],collapse="_"),"\t",tmp.graph[i,2],"\t",tmp.graph[i,1],"\n"),append=TRUE)
    }
    if(verbose){cat("Graph input file for qpGraph written in",outgraphfile,"\n") } 
    ############################
    #######write parameter file for qpGraph
    #######
    parfile=paste0(outfileprefix,".parqpGraph")
    cat(file=parfile,"outpop:  NULL\nforcezmode: YES\nlsqmode: NO\ndiag:  .0001\nbigiter: 6\nhires: YES\nlambdascale: 1\n")
    cat(file=parfile,paste0("fstatsname: ",outfile,"\n"),append=T)
    if(verbose){cat("Parameter File for qpGraph with some default parameters written in",parfile,"\n")}
}