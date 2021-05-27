#' Construct and root an Neighbor-Joining tree of presumably nonadmixed leaves
#' @param fstats Object of class fstats that contains estimates of the fstats (see compute.fstats)
#' @param pop.sel Names of the leaves (pops) used to build the nj tree (at least 3 required)
#' @param edge.fact The multiplying factor of edges length in graph representation 
#' @param plot.nj If TRUE plot the Neighbor-Joining tree
#' @param verbose If TRUE extra information is printed on the terminal 
#' @details A Neighbor-Joining tree is first built (using nj function from the package ape) based on the F2-distance matrix of the leaves in pop.sel which are presumably non-admixed (see the function find.tree.popset to find such groups of scaffold populations using estimated F3 and F4 test statistics). For non-admixed leaves, F2 are indeed expected to be additive along the resulting binary tree (see Lipson et al., 2013). The resulting tree is then rooted using the method described in Lipson et al. (2013) which is based on the property that the estimated heterozygosity of the root h_R equals h_R=1-Q2(A,B) if A and B are two populations sharing R as the only common ancestor in the tree. This estimator should then be consistent across all the possible pairs of populations A and B that are only connected through R in the tree (i.e., that each belong to one of the two partitions of the tree defined by a root position R). Note that 1-Q2(A,B)=(1-Q1(A))/2 + (1-Q1(B))/2 + F2(A,B)=(h_A+h_B)/2+F2(A,B) where h_A, h_B and F2(A,B) are estimated with the function compute.fstats. 
#' @return A list with the following elements:
#' \enumerate{
#' \item "n.rooted.trees": The number of possible rooted binary trees that were evaluated
#' \item "fitted.rooted.trees.list": a list of objects of class fitted.graph containing information on all the possible graphs (indexed from 1 to n.rooted.trees). Each tree may be visualized or further used using functions applied to objects of class fitted.graph (e.g., plot, add.leave) 
#' \item best.rooted.tree The tree (object of class fitted.graph) among all the graphs within fitted.rooted.trees.list displaying the minimal the minimal sd over estimates of h_P (see details) 
#' \item "root.het.est.var": For a matrix of n.tree rows (same order as in the list rooted.tree) and 4 columns with i) the average estimated root heterozygosity h_R across all the pairs of population leave that are relevant for estimation (see details); ii) the size of the range of variation and iii) the s.d. of the estimates of h_R, and iv) the number of population pairs relevant for estimation 
#' \item "nj.tree.eval": If n.edges>3, gives the five worst configuration fit (by calling the compare.fitted.fstats function) which are the same irrespective of rooting
#' }
#' @seealso see \code{\link{fit.graph}}, \code{\link{generate.graph.params}} and \code{\link{add.leaf}}.
#' @export
rooted.njtree.builder<-function(fstats,pop.sel,edge.fact=1000,plot.nj=FALSE,verbose=TRUE){
  if(!is.fstats(fstats)){
    stop("fstats must be an object of class fstats (see the function compute.fstats)\n")
  }else{
    if(!fstats@blockjacknife){
    if(sum(!(pop.sel%in%fstats@comparisons[["pops"]]))>0){
      stop("Some leaves include in pop.sel are not represented in the fstats object\n")
    }
    }
    if(nrow(fstats@heterozygosities)==0){
      stop("fstats object must contain estimates of leaves heterozygosities with standard errors (see computeFstats should be run with compute.heterozygosity=TRUE)\n")  
    }
    if(nrow(fstats@Q.matrix)==0){
        stop("Q.matrix is needed in fstats to compute scores (i.e., compute.fstats must be run with computeF3=TRUE and nsnp.per.bjack.block>0 when generating the fstats object)\n")
        compute.score=FALSE
    }
  }
  n.popsel=length(pop.sel)
  if(n.popsel<3){stop("At least 3 pops are required\n")}
  dum.pops=fstats@comparisons[["pops"]]
  n.pops.fstats.obj=length(dum.pops)
  mat.f2=matrix(0,n.pops.fstats.obj,n.pops.fstats.obj,dimnames=list(dum.pops,dum.pops))
  mat.f2[lower.tri(mat.f2,diag=F)]=fstats@f2.values[,2]
  mat.f2[upper.tri(mat.f2,diag=F)]=t(mat.f2)[upper.tri(mat.f2,diag=F)]
  mat.f2=mat.f2[pop.sel,pop.sel]
  nj.tree=nj(as.dist(mat.f2))
  dist.tree=dist.nodes(nj.tree)
  n.edges=nrow(nj.tree$edge)
  if(plot.nj){plot(nj.tree)}
  if(sum(dist.tree<0)>0){
    stop("Neighbor-Joining tree results in negative branch lengths.\nThe populations are likely not related by a simple binary tree.\n")
  }
  #compute score (for information)
   f2.pred=dist.tree[1:length(nj.tree$tip.label),1:length(nj.tree$tip.label)]
   dimnames(f2.pred)=list(nj.tree$tip.label,nj.tree$tip.label)
   pop.ref=pop.sel[1]
   
   tmp.sel.f2.idx=which((fstats@comparisons[["F2"]][,1]==pop.ref & fstats@comparisons[["F2"]][,2]%in%pop.sel) | (fstats@comparisons[["F2"]][,2]==pop.ref & fstats@comparisons[["F2"]][,1]%in%pop.sel))
   vect.f2.obs=fstats@f2.values[tmp.sel.f2.idx,2]
   tmp.f2.names=fstats@comparisons[["F2"]][tmp.sel.f2.idx,]
   vect.f2.pred=vect.f2.obs*0
   for(i in 1:length(vect.f2.pred)){vect.f2.pred[i]=f2.pred[tmp.f2.names[i,1],tmp.f2.names[i,2]]}
   
   tmp.sel.f3.idx=which(fstats@comparisons[["F3"]][,1]==pop.ref & fstats@comparisons[["F3"]][,2]%in%pop.sel  & fstats@comparisons[["F3"]][,3]%in%pop.sel)
   vect.f3.obs=fstats@f3.values[tmp.sel.f3.idx,2]
   if(length(tmp.sel.f3.idx)==1){#evite pb si n=1
     tmp.f3.names=matrix(fstats@comparisons[["F3"]][tmp.sel.f3.idx,],ncol=3)
     rownames(tmp.f3.names)=rownames(fstats@comparisons[["F3"]])[tmp.sel.f3.idx]
   }else{
    tmp.f3.names=fstats@comparisons[["F3"]][tmp.sel.f3.idx,]
   }
   vect.f3.pred=vect.f3.obs*0
   for(i in 1:length(vect.f3.pred)){#F3(A;B,C)=(F2(A,B)+F2(A,C)-F2(B,C))/2
     vect.f3.pred[i]=(f2.pred[tmp.f3.names[i,1],tmp.f3.names[i,2]]+f2.pred[tmp.f3.names[i,1],tmp.f3.names[i,3]]-f2.pred[tmp.f3.names[i,2],tmp.f3.names[i,3]])/2
   }
   tmp.names=c(rownames(tmp.f2.names),rownames(tmp.f3.names))
   tmp.Qmat=fstats@Q.matrix[tmp.names,tmp.names]
   tmp=c(vect.f2.pred,vect.f3.pred)-c(vect.f2.obs,vect.f3.obs)
  # ls.score=sum(tmp)**2
  # cat("LS score of the NJ tree:",ls.score,"\n")
   wls.score=as.numeric(t(tmp)%*%solve(tmp.Qmat)%*%tmp)
   if(verbose){cat("Score of the NJ tree:", wls.score,"\n")} 
   #bic for informtaion
   Kvalue=length(tmp)*log(2*pi)+determinant(tmp.Qmat,logarithm=TRUE)$modulus[1] #-2*cte de likelihood
   bic=wls.score+n.edges*log(length(tmp)) - Kvalue
   fitted.outstats=cbind(c(vect.f2.pred,vect.f3.pred),c(vect.f2.obs,vect.f3.obs),tmp/sqrt(diag(tmp.Qmat)))
   colnames(fitted.outstats)=c("Stat. value","Fitted Value","Z-score")
   row.names(fitted.outstats)=tmp.names
  #computation of hP (heterozygosities at ancestral populaion of each pairs of pops: cf. Lipson criterion to root trees)
  #hP=1-Q2(A,B)=1-[(Q1(A)+Q1(B))/2 - F2(A,B)]=(1-Q1(A))/2 + (1-Q1(B))/2 + F2(A,B) #avoid to store Q2 in fstats
  dum=mat.f2*0 + fstats@heterozygosities[pop.sel,2]
  hP=mat.f2 + (dum+t(dum))/2
  ###
  rooted.graph=list() # for each possible position of the root give the resulting graphs in appropriate format
  rooted.graph.edges.length=list() #for each graph the root is put in the middle of the candidate tip
  if(verbose){
  cat("Construction of all the",n.edges,"possible rooted tree from the NJ tree\n(stored as graph in the rooted.graph object of the output list)\n")
  }
  n.cdt.edges=nrow(nj.tree$edge)
  n.tips=length(nj.tree$tip.label)
  tip.list=1:n.tips
  nroot=nj.tree$Nnode + n.tips + 1
  for(i in 1:n.edges){
    posroot=i
    edgeinfo=nj.tree$edge[posroot,]
    #First create the two edges from the root
    newedges=rbind(c(nroot,edgeinfo[1],nj.tree$edge.length[posroot]/2),c(nroot,edgeinfo[2],nj.tree$edge.length[posroot]/2))
    #go down from root to tips
    oldedges=cbind(nj.tree$edge[-posroot,],nj.tree$edge.length[-posroot])
    parents=newedges[,2]
    while(length(parents)>0){
      parents.new=c() #initialize to allow exiting the while loop
      # to the left
      child=which(oldedges[,1] %in% parents)
      if(length(child)>0){
        parents.new=oldedges[child,2]
        newedges=rbind(newedges,oldedges[child,1:3])
        oldedges=oldedges[-child,]
      }
      #a droite
      child=which(oldedges[,2] %in% parents)
      if(length(child)>0){
        parents.new=c(parents.new,oldedges[child,1])
        newedges=rbind(newedges,oldedges[child,c(2:1,3)])
        oldedges=oldedges[-child,]
      }
      parents=parents.new
    }
    node.names=c(nj.tree$tip.label,paste0("I",1:nj.tree$Nnode),"Root")
    rooted.graph[[i]]=cbind(node.names[newedges[,2]],node.names[newedges[,1]],rep("",nrow(newedges)))
    rooted.graph.edges.length[[i]]=newedges[,3]
  }
  
  ###
  ###Compare root position based on the variation in root heterozygosity estimates among all pairs of populations from each side of th root
  comp.hp=matrix(0,n.edges,4)
  colnames(comp.hp)=c("Mean","Range","sd","ncomps")
  for(g in 1:n.edges){
    tmp.graph=rooted.graph[[g]] #by construction the two root edges coming from the root are the first ones
    #left partition
    tmp.left.part=tmp.graph[1,1]
    cur.par=tmp.left.part ; crit=T
    while(crit){
      cur.edge=tmp.graph[,2]%in%cur.par
      if(sum(cur.edge)>0){
        cur.par=tmp.graph[cur.edge,1]
        tmp.left.part=c(tmp.left.part,cur.par)
      }else{crit=F}
    }
    tmp.left.part=tmp.left.part[tmp.left.part%in%pop.sel]
    #right partition
    tmp.right.part=tmp.graph[2,1]
    cur.par=tmp.right.part ; crit=T
    while(crit){
      cur.edge=tmp.graph[,2]%in%cur.par
      if(sum(cur.edge)>0){
        cur.par=tmp.graph[cur.edge,1]
        tmp.right.part=c(tmp.right.part,cur.par)
      }else{crit=F}   
    }
    tmp.right.part=tmp.right.part[tmp.right.part%in%pop.sel]   
    
    dd=hP[tmp.left.part,tmp.right.part]
    comp.hp[g,]=c(mean(dd),diff(range(dd)),sd(dd),length(dd))
  }
  
  ### construction de l'objet fitted graph 
  if(nrow(fstats@heterozygosities)>0){
    drift.scaling=TRUE
    hetero=fstats@heterozygosities[pop.sel,2]
  }
  for(g in 1:n.edges){
    dum=new("fitted.graph")
    tmp.graph.params=generate.graph.params(rooted.graph[[g]])
    edges.length=rooted.graph.edges.length[[g]]
    names(edges.length)=colnames(tmp.graph.params@graph.matrix)
    dum@graph=tmp.graph.params@graph
    dum@score=wls.score
    dum@bic=bic
    dum@fitted.outstats=fitted.outstats
    dum@edges.length=edges.length
    #computed Omega
    suppressWarnings(graph.path.to.root<-matrix(as.numeric(tmp.graph.params@graph.matrix),nrow = tmp.graph.params@n.leaves,ncol=ncol(tmp.graph.params@graph.matrix) ,dimnames = dimnames(tmp.graph.params@graph.matrix)))
    tmp.omega=t(t(graph.path.to.root)*dum@edges.length)
    tmp.omega=tmp.omega%*%t(graph.path.to.root)
  #  dum@fitted.omega=tmp.omega
    f2.mat=diag(tmp.omega)%*%t(rep(1,nrow(tmp.omega)))
    f2.mat=f2.mat+t(f2.mat)-2*tmp.omega
    dum@fitted.f2.mat=f2.mat
    ###Drift scaling?
    if(drift.scaling){
      tmp.edges=matrix(unlist(strsplit(names(edges.length),split="<->")),nrow=length(edges.length),byrow=T)
      tmp.nodes=unique(as.vector(tmp.edges))
      tmp.het=rep(NA,length(tmp.nodes))
      names(tmp.het)=tmp.nodes
      tmp.het[pop.sel]=hetero
      for(i in pop.sel){
        if(i %in% tmp.edges[,2]){
         crit=T ; tmp.cur=i
          while(crit){
              tmp.e=tmp.edges[,2]==tmp.cur
              if(sum(tmp.e)==0){
                crit=F
              }else{
                tmp.p=tmp.edges[tmp.e,1]
                tmp.het[tmp.p]=tmp.het[i]+2*edges.length[tmp.e]
                tmp.cur=tmp.p
              }
            }
          }}
        if(sum(is.na(tmp.het))>0){
          if(verbose){
           cat("Warning: Problem in estimating internal node heterozygosity\n")
           cat("No scaling will be done\n")
          }
          drift.scaling=FALSE
        }else{
          edges.length.scaled=2*edges.length/tmp.het[tmp.edges[,1]]
          nodes.het=tmp.het
          dum@edges.length.scaled=edges.length.scaled
        }
    }
    ####graph representation of the results
    ##header of the digraph file
    leaves.color="green";size="7.5,10" #for dot file
    outlines=c("digraph G {\n",paste0("size = \"",size,"\" ;\n"),"labelloc = \"t\" ;\n","")
    ##info leaves
    leaves.color=rep(leaves.color,length(pop.sel))
    outlines=c(outlines,paste0("\"",pop.sel,"\" [ label = \"",pop.sel,"\",color=",leaves.color,",style=filled ] ;\n"),"")
    ##info edges
    tmp.nodes=matrix(unlist(strsplit(names(edges.length),split="<->")),ncol=2,byrow = T)
    if(drift.scaling){
      tmp.values=as.character(round(edges.length.scaled*edge.fact))
    }else{
      tmp.values=as.character(round(edges.length*edge.fact))
    }
    tmp.info=rep("style=plain",nrow(tmp.nodes))
    tmp=paste0("\"",tmp.nodes[,1],"\"->\"",tmp.nodes[,2],"\" [ label=\"",tmp.values,"\",",tmp.info," ] ;\n")
    outlines=c(outlines,tmp,"}")
    dum@dot.graph=outlines

    rooted.graph[[g]]=dum
              }

  ###
  best.tree=rooted.graph[[which.min(comp.hp[,3])]]
  rownames(comp.hp)=paste0("Tree",1:nrow(comp.hp))
  ###eval fit (worst Z-score)
  if(n.edges>3){#quels que soit les graphe on obtient la meme chose! (on prend le best.tree en ref)
    tmp.comp=compare.fitted.fstats(fstats,best.tree,n.worst.stats=0)
    tmp.which.worst=order(abs(tmp.comp[,3]),decreasing = T)[1:5]
    eval.det=tmp.comp[tmp.which.worst,]
    return(list(n.rooted.trees=n.edges,fitted.rooted.trees.list=rooted.graph,best.rooted.tree=best.tree,root.het.est.var=comp.hp,nj.tree.eval=eval.det))    
  }else{
   return(list(n.rooted.trees=n.edges,fitted.rooted.trees.list=rooted.graph,best.rooted.tree=best.tree,root.het.est.var=comp.hp))
  }
}
