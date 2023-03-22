#' Generate a graph parameter object to fit admixture graph to observed fstats
#' @param graph A three columns matrix containing graph information in a simple format (see details)
#' @param fstats A fstats object containing estimates of fstats
#' @param popref Reference population of the fstats basis used to fit the graph.
#' @param outfileprefix The prefix of the dot file that will represent the graph (with extension ".dot"). If NULL, no graph file generated
#' @param verbose If TRUE some information is printed on the terminal
#' @details The graph needs to be specified by a three column (character) matrix corresponding for each edge (wether admixed or not) to i) the child node; ii) the parent node; iii) the admixture proportion. For non-admixed edge, the third column must be blank. An admixed node should be referred two times as a child node with two different parent node and two different admixture proportions coded as alpha and (1-alpha) (Note that the parentheses are mandatory) if alpha is the name of the admixture proportion. The root is automatically identified as a node only present in the parent node column. Several checks are made within the function but it is recommended to check the graph by plotting the resulting dot file named {outfileprefix.}dot using for instance the grViz() from the DiagrammeR package that may be called directly with plot or with the dot program (e.g., dot -Tpng inputgraph.dot in terminal). Note that the dot file may be easily customized (e.g., to change leave color, parameter names...). The fstats object should be of class fstats (see help(fstats) for details) containing estimates of F2 and F3 statistics and block jackknife as generated with the \code{\link{compute.fstats}} function with computeF3 set to TRUE. If no fstats object is provided, only graph parameters will be generated.  
#' @return An object of class graph.params (see help(graph.params) for details)
#' @seealso The object may be used to estimate graph parameters with the function \code{\link{fit.graph}} or to generate files for the qpGraph software with \code{\link{graph.params2qpGraphFiles}}. See also \code{\link{graph.params2symbolic.fstats}} to obtain symbolic representation of Fstats.
#' @examples
#' graph=rbind(c("P1","P7",""),c("P2","s1",""),c("P3","s2",""),c("P6","S",""),
#'             c("S","s1","a"),c("S","s2","(1-a)"),c("s2","P8",""),c("s1","P7",""),
#'             c("P4","P9",""),c("P5","P9",""),c("P7","P8",""),
#'             c("P8","R",""),c("P9","R","")) 
#' graph.params=generate.graph.params(graph)
#' plot(graph.params)
#' ##NOTE: this calls grViz from DiagrammeR which cannot easily be plotted 
#' #within pdf or other device. To that end the easiest is to output 
#' #the graph in a dot file (using the outfileprefix argument) and 
#' #then to use the dot program out of R in a terminal: dot -Tpng inputgraph.dot
#' @export
generate.graph.params<-function(graph,fstats=NULL,popref=NULL,outfileprefix=NULL,verbose=TRUE){
  ##check if any pop or adm parm is named I or i (otherwise may be interpreted as complex number in literal computation)
  if(sum(graph=="I"|graph=="i")>0){
    stop("Pops or Admixture proportion parameters cannot be named I or i (you may change into any another letter(s))\n")
  }
  ##
  out=new("graph.params")
  #######################
  ####prepare graph file in dot format (allow evaluating input graph)
  ######################
  leaves=unique(graph[,1][!(graph[,1]%in%graph[,2])]) #in case a leave is admixed
  n.leaves=length(leaves)
  if(n.leaves<2){stop("Less than two leaves detected (You may check the graph object)\n")}
  root=unique(graph[,2][!(graph[,2]%in%graph[,1])])
  if(length(root)>1){
    cat("More than one root detected in the graph with names:\n")
    show(root)
    stop("The graph should have one and only one root\n")
  }
  if(length(root)==0){stop("No root detected in the graph\n")}
  adm.graph.rows=which(nchar(graph[,3])>0)
  if(length(adm.graph.rows)>0){adm.graph=TRUE}else{adm.graph=FALSE}
  if(adm.graph){
    adm.pops=table(graph[adm.graph.rows,1])
    n.adm.nodes=length(adm.pops)
    if(sum(adm.pops!=2)>0){
      cat("Admixed Pop detected (names and number of parents):")
      show(adm.pops)
      stop("All admixed pop should have two parental pop\n")
    }
    #check noms params
    dum=graph[adm.graph.rows,3]
    dum.eval=!grepl("^[a-zA-Z0-9()-]*$",dum)
    if(sum(dum.eval)>0){
      cat("PB: forbidden characters detected in admixture proportions:\n")
      show(dum[dum.eval])
      stop("Only alphanumeric characters tolerated in admixture proportion names\n")
    }
    #check si plusieurs signes moins (au cas ou nom admixt params avec un signe -)
    dum.eval=(nchar(dum)-nchar(gsub("-","",dum)))>1
    if(sum(dum.eval)>0){
      cat("PB: More than one minus sign - detected in an admixture proportion:\n")
      show(dum[dum.eval])
      stop("Only alphanumeric characters tolerated in admixture names\n")
    }
    #match des parentheses: on met une erreur car trop complique a gerer sinon
    dum.eval=(nchar(gsub("[(]","",dum))-nchar(gsub("[)]","",dum)))!=0 |
             (nchar(gsub("[(]","",dum))-nchar(gsub("-","",dum)))!=0 |
             (nchar(gsub("[(]","",dum))-nchar(gsub("-","",dum)))!=0 
    if(sum(dum.eval)>0){
      cat("PB: check matching (or presence) of parentheses in:\n")
      show(dum[dum.eval])
      stop("admixture proportion must be of the form alpha and (1-alpha) and parentheses are mandatory in the latter case\n")
    }
  }else{
    n.adm.nodes=0
  }
  ####################
  ##header of the digraph file
  leaves.color="green";size="7.5,10" #for dot file
  outlines=c("digraph G {\n",paste0("size = \"",size,"\" ;\n"),"labelloc = \"t\" ;\n","")
  ##info leaves
  leaves.color=rep(leaves.color,n.leaves)
  outlines=c(outlines,paste0("\"",leaves,"\" [ label = \"",
                             leaves,"\",color=",leaves.color,",style=filled ] ;\n"),"")
  ##info edges (non admixed)
  ##admixed edges are set to null (as used to solve identifiability issue: cf nullify)
  n.edges=nrow(graph)-length(adm.graph.rows)
  tmp.info=rep("style=plain",nrow(graph)) ;  tmp.info[adm.graph.rows]="style=dotted"
  tmp=paste0("\"",graph[,2],"\"->\"",graph[,1],"\" [ label=\"",graph[,3],"\",",tmp.info," ] ;\n")
  outlines=c(outlines,tmp,"}")
  out@dot.graph=outlines
  if(!is.null(outfileprefix)){
   writeLines(outlines,con=paste0(outfileprefix,".inputgraph.dot"))
   cat("Graph input file in dot format written in",paste0(outfileprefix,".inputgraph.dot"),".\nIt can be visualized using the grViz function from the DiagrammeR package, i.e., by running the command:",paste0("grViz(\"",outfileprefix,".inputgraph.dot\")"),"\n")
  }
  #####################
  ##Compute paths_to_root_matrix
  #####################
  nodes=unique(as.vector(graph[,1:2])) ; n.nodes=length(nodes) #all nodes (leaves+internal nodes + root)
  root_idx=which(nodes==root)
  edges=graph[,2:1] #parent, descendant
  n.edges=nrow(graph) #admixture edges are also counted will need to be nullified afterwards (cf. à la treemix ou mixmapper)
  edges.names=paste(graph[,2],graph[,1],sep="<->")
  rownames(edges)=edges.names
  edges_idx=1:n.edges
  edges_idx.matrix=matrix(NA,n.nodes,n.nodes,dimnames=list(nodes,nodes)) #pour gere la direction des chemin (from to)
  for(i in 1:n.edges){edges_idx.matrix[graph[i,1],graph[i,2]]=edges_idx.matrix[graph[i,2],graph[i,1]]=i}
  probs=matrix("",n.nodes,n.nodes,dimnames=list(nodes,nodes)) #proba admixture
  parents=matrix(FALSE,n.nodes,n.nodes,dimnames=list(nodes,nodes)) 
  for(i in 1:nrow(graph)){
    probs[graph[i,1],graph[i,2]]=probs[graph[i,2],graph[i,1]]=graph[i,3]
    parents[graph[i,1],graph[i,2]]=TRUE
  }
  paths_to_root_matrix=matrix("",n.leaves,n.edges,dimnames = list(leaves,edges.names))
  root.edges=graph[,2]==root #useful to find idx of edges in the final paths_to_root_matrix (after nullifying admixture edges if any)
  for(src in leaves){
    npaths=1 ; proba.paths="1"
    src_idx <- which(src == rownames(parents))    
    cur.node=src_idx
    while(npaths>0){
      cur.node.new=proba.paths.new=c()
      for(i in 1:npaths){
        parent.node=which(parents[cur.node[i],])
        cur.node.new=c(cur.node.new,parent.node)
        for(j in 1:length(parent.node)){
          tmp.prob=proba.paths[i]
          if(nchar(probs[cur.node[i],parent.node[j]])>0){
            tmp.prob=probs[cur.node[i],parent.node[j]]
            if(proba.paths[i]!="1"){tmp.prob=paste(proba.paths[i],tmp.prob,sep="*")}
          }
          proba.paths.new=c(proba.paths.new,tmp.prob)
          tmp.edge.id=edges_idx.matrix[cur.node[i],parent.node[j]]
          tmp.prob.mat=paths_to_root_matrix[src,tmp.edge.id]
          if(nchar(tmp.prob.mat)>0){tmp.prob.mat=paste(tmp.prob.mat,tmp.prob,sep="+")}else{tmp.prob.mat=tmp.prob}
          paths_to_root_matrix[src,tmp.edge.id]=tmp.prob.mat
        }
      }
      tmp.eval=cur.node.new!=root_idx
      npaths=sum(tmp.eval)
      if(npaths>0){
        cur.node=cur.node.new[tmp.eval]
        proba.paths=proba.paths.new[tmp.eval]
      }
    }
  }
  paths_to_root_matrix[paths_to_root_matrix==""]="0"
  ## Mise à 0 des edges menant à une pop admixee (gestion des pb de non-identifiabilite des parametre, cf treemix et Lipson) et fonction nullify_admixtpop_source_edges
  ## Exemple: O1---S1-(alpha)---A--(1-alpha)--S2----O2
  ##                            |
  ##                            T
  ## On considere O1--S1=O2--S2=0 (Rk: dans TreeMix, c'est A-T=0 et O1--S1=0 (si alpha<0.5) ou O2--S2=0 (si alpha>0.5))
  # cf   nullify_admixtpop_source_edges<-function(graph,graph.path.to.root)
  edges.admixed=nchar(graph[,3])>0  
  paths_to_root_matrix=paths_to_root_matrix[,!edges.admixed]
  ##simplification des expressions avec yacas
  for(i in 1:nrow(paths_to_root_matrix)){
    for(j in 1:ncol(paths_to_root_matrix)){
      if(nchar(paths_to_root_matrix[i,j])>1){
        paths_to_root_matrix[i,j]=Ryacas::yac(paste0("Expand(",paths_to_root_matrix[i,j],")"),rettype="str")
      }
    }
  }

#remplissage de l'objet  
  out@edges.names=colnames(paths_to_root_matrix)
  out@n.edges=ncol(paths_to_root_matrix)  
  out@graph.matrix=paths_to_root_matrix
  out@leaves=leaves ;  out@n.leaves=n.leaves
  out@n.adm.nodes=n.adm.nodes  ;  out@graph=graph
  out@is.admgraph=adm.graph ;  out@n.nodes=n.nodes
  out@nodes.names=nodes ;  out@root.name=root  
  root.edges=root.edges[!edges.admixed]
  out@root.edges.idx=which(root.edges)  
  if(adm.graph){
    tmp=unique(gsub("[-1()]","",graph[edges.admixed,3]))
    if(length(tmp)!=n.adm.nodes){stop("Problem in identifying admixture parameter names : should be of the form alpha and (1-alpha): parentheses are important in the latter case\n")}
    out@adm.params.names=tmp 
  }
  
  #######################
  ####retrieve F2 and F3 stats for the leaves: relevant F2 and F3: (nleaves-1) F2 (popref with all tohers) + (nleaves-2)*(nleaves-3)/2 F3 involving popref as target
  ####Pour construire x.mat symbolique: si fstats on utilise ordre des comparisons F2 et F3 definis dans objets fstats (qui depend de popnames objet data alors que f2.symbolic depend de leaves defini dans graph); si pas fstats on utilise ordre des leaves (qui depend de graph)=> on peut obtenir un ordre des lignes de x.mat selon qu'on donne fstats ou pas masi c'est normal!
  ######################
  if(!is.null(fstats)){
    ##recuperation popref
    if(is.null(popref)){#popref may be useful to compute equation (was only defined before if fstats non null)
      popref=leaves[1]
    }else{
      if(!(popref %in% leaves)){
        warning("popref was not found among the leaves (the leave will be used instead")
        popref=leaves[1]
      }
    }
    out@popref=popref
    ##some checks
    if(!is.fstats(fstats)){
      stop("fstats must be an object of class fstats (see the function compute.fstats)\n")
    }
    if(!fstats@blockjacknife){
        stop("No block-jackknife estimate of estimators s.e. found in the object fstats.\nFunctions compute.fstats must be run with nsnp.per.block>0\n")
    }
    if(nrow(fstats@f3.values)==0){
      stop("No estimates of F3 statistics available in the fstats object.\nFunctions compute.fstats must be run with computeF3=TRUE\n")
    }
    #sel F2 stats
    n.f2.exp=n.leaves*(n.leaves-1)/2
    dum.sel=which((fstats@comparisons[["F2"]][,1]%in%leaves) & (fstats@comparisons[["F2"]][,2]%in%leaves))
    if(length(dum.sel)!=n.f2.exp){stop("Some combination of F2 for leaves are not present in fstats$F2\n")}
    f2.idx=dum.sel[rowSums(fstats@comparisons[["F2"]][dum.sel,1:2]==popref)==1]
    f2.pops=fstats@comparisons[["F2"]][f2.idx,]
    ##for practical reason: check f2.pop to always have refpop in first position (Rk: F2 is symetric=> no change of value)
    tmp.mod=f2.pops[,1]!=popref
    if(sum(tmp.mod)>0){f2.pops[tmp.mod,2]=f2.pops[tmp.mod,1] ; f2.pops[,1]=popref}
    out@f2.target=fstats@f2.values[f2.idx,2]
    rownames(f2.pops)=paste0(f2.pops[,1],",",f2.pops[,2]) #replace P1 in first position for names (used in fit.graph to fill mat.f2)
    out@f2.target.pops=f2.pops
    #sel F3 stats    
    n.f3.exp=(n.leaves-1)*(n.leaves-2)/2 #i.e., those involving popref as target
    dum.pops=fstats@comparisons[["F3"]]
    f3.idx=which(dum.pops[,1]==popref & (dum.pops[,2]%in%leaves) & (dum.pops[,3]%in%leaves))
    if(length(f3.idx)!=n.f3.exp){stop("Some combination of F3 for leaves are not present in fstats$F3\n")}
    f3.pops=dum.pops[f3.idx,]
    out@f3.target=fstats@f3.values[f3.idx,2]
    out@f3.target.pops=f3.pops
    ######exporting f2, f3, Qmat...
    tmp.sel.qmat=c(rownames(fstats@f2.values)[f2.idx],rownames(fstats@f3.values)[f3.idx])
    out@f.Qmat=fstats@Q.matrix[tmp.sel.qmat,tmp.sel.qmat]
    ##############################
    if(verbose){
      N.edges.length=ncol(out@graph.matrix)-1 #unrooted
      tmp.nparams=N.edges.length+n.adm.nodes
      tmp.nstats=nrow(f2.pops)+nrow(f3.pops)
      cat("Total Number of Parameters:",N.edges.length+n.adm.nodes,paste0("(",N.edges.length," edges lengths + ",n.adm.nodes," adm. coeff.)"),"\n") #-1 because the two root branches are not identifiable (unrooted if using F2...)
     cat("Total Number of Statistics:",tmp.nstats,paste0("(",nrow(f2.pops)," F2 and ",nrow(f3.pops)," F3)"),"\n")
     if(tmp.nparams>tmp.nstats){
       cat("WARNING: Overfitting: the number of parameters exceeds the number of statistics. The fit.graph function will return a null value on the obtained graph.params object.\n")
     }
    }
    ######retrieve heterozygosity if available
    if(nrow(fstats@heterozygosities)>0){#forcement jackknif si on arrive jusqu'ici
      het=fstats@heterozygosities[,2] 
      names(het)=rownames(fstats@heterozygosities)
      het=het[leaves]
      out@Het=het
    }
  }
  ###########
  return(out)
}
