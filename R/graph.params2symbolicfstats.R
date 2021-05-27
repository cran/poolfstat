#' Provide a symbolic representation of all the F-statistics and the model system of equations
#' @param x An object of class graph.params containing graph information and relevant Fstats estimates (see the function generate.graph.params)
#' @param outfile The file where to print the equations (default=NULL, equations are not printed in a file)
#' @return A list with the following elements:
#' \enumerate{
#' \item "model.matrix": A symbolic representation of the matrix M relating the basis F-statistics and graph edge length as F=M*b where F is the vector of the basis Fstats (row names of model.matrix M) and b is the vector of graph edges (column names of model.matrix M).
#' \item "omega": A symbolic representation of the scaled covariance matrix of allele frequency with edge names and admixture parameter names as specified in the edges.names and adm.params.names slot of the input graph.params object x
#' \item "F2.equations": A symbolic representation of the nleaves(nleaves-1)/2 different F2 as a function of graph parameters
#' \item "F3.equations": A symbolic representation  of the nleaves(nleaves-1)(nleaves-2)/2 different F3 as a function of graph parameters 
#' \item "F4.equations": A symbolic representation of the npops(npops-1)(npops-2)(npops-3)/8 different F4 as a function of graph parameters
#' }
#' @seealso To generate a graph.params object, see \code{\link{generate.graph.params}}. 
#' @examples
#' graph=rbind(c("P1","P7",""),c("P2","s1",""),c("P3","s2",""),c("P6","S",""),
#'             c("S","s1","a"),c("S","s2","(1-a)"),c("s2","P8",""),c("s1","P7",""),
#'             c("P4","P9",""),c("P5","P9",""),c("P7","P8",""),
#'             c("P8","R",""),c("P9","R","")) 
#' graph.params=generate.graph.params(graph)
#' graph.equations=graph.params2symbolic.fstats(graph.params) 
#' @export

graph.params2symbolic.fstats<-function(x,outfile=NULL){
  if(!is.graph.params(x)){stop("The input x is not a valid graph.params object (see the function generate.graph.params)\n")}
  tmp.edges.code=paste0("zZz",1:x@n.edges,"zZz") #on recode pour empecher gag dans calcul symbolique (par exemple des caracteres comme _ @ posent probleme)
  if(!is.character(outfile)){
    printfile=FALSE
  }else{
    printfile=TRUE
    cat("Equations will be printed in file",outfile,"\n")
  }  
  ########################
  ###Symbolic representation of Omega
  ########################
  omega=matrix("0",x@n.leaves,x@n.leaves,dimnames=list(x@leaves,x@leaves))
  paths.with.parenthesis=x@graph.matrix
  #il faut mettre des () dans les cas ou on additionne des poids (e.g., a+(1-a) ou a*b + (1-a)*b...: dans ces on a plus de signes + ou - que de parentheses ouvrantes)
  tmp.add.parenthesis=nchar(gsub("[-+]","",x@graph.matrix))<nchar(gsub("\\(","",x@graph.matrix)) #les "-" des coef admixt annule les ( des coef admixt)
  paths.with.parenthesis[tmp.add.parenthesis]=paste0("(",paths.with.parenthesis[tmp.add.parenthesis],")")
  for(i in 1:x@n.leaves){
    tmp.sel1=x@graph.matrix[i,]!="0"
    for(j in i:x@n.leaves){
      tmp.sel2=x@graph.matrix[j,]!="0"
      tmp.sel=tmp.sel1 & tmp.sel2
      if(sum(tmp.sel)>0){
        tmp.val=paste(paths.with.parenthesis[i,tmp.sel],paths.with.parenthesis[j,tmp.sel],tmp.edges.code[tmp.sel],sep="*")
        ##simplification avec yacas
        tmp.val=Ryacas::yac(paste0("Expand(",paste(tmp.val,collapse="+"),")"),rettype="str")
        omega[i,j]=omega[j,i]=tmp.val
      }
    }
  }
  #on recode les noms
  for(i in 1:x@n.edges){omega=gsub(tmp.edges.code[i],x@edges.names[i],omega)}

  #######################
  ##Construction F2 (edge names recodes)
  #######################
  tmp.omega=omega
  for(i in 1:x@n.edges){tmp.omega=gsub(x@edges.names[i],tmp.edges.code[i],tmp.omega)}
  mat.f2=matrix("0",x@n.leaves,x@n.leaves)
  for(i in 1:(x@n.leaves-1)){
    for(j in (i+1):x@n.leaves){
      tmp=paste0(tmp.omega[i,i]," + ",tmp.omega[j,j]," - 2*(",tmp.omega[i,j],")")
      mat.f2[i,j]=mat.f2[j,i]=Ryacas::yac(paste0("Expand(",tmp,")"),rettype="str")
    }
  }
  dimnames(mat.f2)=dimnames(omega)
  
  #######################
  ##Construction de x.mat et impression
  #######################  
  if(length(x@f2.target)==0 | length(x@f3.target)==0){
      #calcul de x.mat avec ordre correspond a leaves (pas besoin d etre sur de la correspondance avec les stats)
      popref=x@leaves[1]
      popref.idx=which(rownames(mat.f2)==popref)
      vect.fstats=mat.f2[popref.idx,][-popref.idx]
      vect.fstats.names=paste0("F2(",popref,",",rownames(mat.f2)[-popref.idx],")")
      alt.pop.idx=(1:nrow(mat.f2))[-popref.idx]
      for(i in 1:(length(alt.pop.idx)-1)){
        for(j in (i+1):length(alt.pop.idx)){
          tmp=paste0("(",mat.f2[popref.idx,alt.pop.idx[i]],"+",mat.f2[popref.idx,alt.pop.idx[j]],"-(",mat.f2[alt.pop.idx[i],alt.pop.idx[j]],"))/2")
          tmp=Ryacas::yac(paste0("Expand(",tmp,")"),rettype="str")
          vect.fstats=c(vect.fstats,tmp)
          vect.fstats.names=c(vect.fstats.names,paste0("F3(",popref,";",rownames(mat.f2)[alt.pop.idx[i]],",",rownames(mat.f2)[alt.pop.idx[j]],")"))
        }
      }
   }else{
  ##Si objet fstat disponible: on respecte bien les meme stats frouni dans objet fstats,; cf infra si pas objet fstats)
  ##Remise dans le meme ordre de x.mat (a verifier car pas la meme procedure de recuperation des F2)
    popref=x@popref
    leaves=rownames(mat.f2)
    popref.idx=which(leaves==popref)
    vect.fstats=c()
    #F2
    for(i in 1:nrow(x@f2.target.pops)){vect.fstats=c(vect.fstats,mat.f2[popref.idx,which(leaves==x@f2.target.pops[i,2])])}
    #F3
    for(i in 1:nrow(x@f3.target.pops)){
     dum.idx1=which(leaves==x@f3.target.pops[i,2])
     dum.idx2=which(leaves==x@f3.target.pops[i,3])  
     tmp=paste0("(",mat.f2[popref.idx,dum.idx1],"+",mat.f2[popref.idx,dum.idx2],"-(",mat.f2[dum.idx1,dum.idx2],"))/2")
    tmp=Ryacas::yac(paste0("Expand(",tmp,")"),rettype="str")
    vect.fstats=c(vect.fstats,tmp)
  }
  vect.fstats.names=c(paste0("F2(",rownames(x@f2.target.pops),")"),paste0("F3(",rownames(x@f3.target.pops),")"))   
  }
  x.mat=matrix("",length(vect.fstats),x@n.edges)
  for(i in 1:x@n.edges){
    for(j in 1:length(vect.fstats)){
      x.mat[j,i]=Ryacas::yac(paste0("Coef(",vect.fstats[j],",",tmp.edges.code[i],",1)"))
    } 
  }
  dimnames(x.mat)=list(vect.fstats.names,x@edges.names) 

  ######################
  ###system equations
  #####################
  ##creation du systeme dequation a imprimer
  dd=nrow(x.mat) ; pos.signs=floor(dd/2)
  v.egal=v.mult=rep("",dd) ;  v.egal[pos.signs]="=" ;  v.mult[pos.signs]="X"
  tmp.parenthesis=rep("|",nrow(x.mat))
  obj.to.print=cbind(tmp.parenthesis,rownames(x.mat),tmp.parenthesis,v.egal,tmp.parenthesis,x.mat,tmp.parenthesis,v.mult)
  if(nrow(x.mat)>x@n.edges){
    tmp.ne=floor((nrow(x.mat)-x@n.edges)/2) 
    tmp.nb=nrow(x.mat)-x@n.edges-tmp.ne 
    tmp.parenthesis=c(rep("",tmp.nb),tmp.parenthesis[1:x@n.edges],rep("",tmp.ne))
    tmp.val=c(rep("",tmp.nb),x@edges.names,rep("",tmp.ne))
    obj.to.print=cbind(obj.to.print,tmp.parenthesis,tmp.val,tmp.parenthesis)
  }else{
    if(nrow(x.mat)==x@n.edges){
      obj.to.print=cbind(obj.to.print,tmp.parenthesis,x@edges.names,tmp.parenthesis) 
    }else{#ca devrait jamais arriver car overparametre!
      tmp.ne=floor((x@n.edges-nrow(x.mat))/2)
      tmp.nb=x@n.edges-nrow(x.mat)-tmp.ne
      obj.to.print=rbind(
        matrix("",tmp.nb,ncol(obj.to.print)),
        obj.to.print,
        matrix("",tmp.nb,ncol(obj.to.print))  )
      tmp.parenthesis=c(rep("|",tmp.nb),tmp.parenthesis,rep("|",tmp.ne))    
      obj.to.print=cbind(obj.to.print,tmp.parenthesis,x@edges.names,tmp.parenthesis)  
    } 
  }
  ##preparation du vecteur a imprimer (en writeLines pour centrer...)
  tmp.nchar=nchar(obj.to.print)
  for(i in 1:ncol(obj.to.print)){
    tmp.max=max(tmp.nchar[,i])
    tmp.diff=tmp.max-tmp.nchar
    tmp.add.left=floor((tmp.max-tmp.nchar[,i])/2)
    tmp.add.right=tmp.max-tmp.nchar[,i]-tmp.add.left
    tmp.new=cbind(
      unlist(mapply(f<-function(c){paste0(rep(" ",c),collapse = "")},tmp.add.left)),
      obj.to.print[,i],
      unlist(mapply(f<-function(c){paste0(rep(" ",c),collapse = "")},tmp.add.right)))
    obj.to.print[,i]=apply(tmp.new,1,paste,collapse="")
  }
  vect.to.print=apply(obj.to.print,1,paste,collapse=" ")

  if(printfile){  
  cat(file=outfile,"#########################################\n########### Model Equations #############\n#########################################\n\n")
  write(vect.to.print,file=outfile,append = TRUE)
  }
  #################
  #####F3
  ################  
  all.F3.eqs=c() #brute force since character (and fast enough (don't expect lots of pops!))
  for(i in 1:x@n.leaves){
    for(j in 1:(x@n.leaves-1)){
      for(k in (j+1):x@n.leaves){
        if(i!=j & i!=k){
          tmp.nom=paste0("F3(",x@leaves[i],";",x@leaves[j],",",x@leaves[k],") = ")
          tmp.val=paste0("(",mat.f2[i,j],"+",mat.f2[i,k],"-(",mat.f2[j,k],"))/2")
          tmp.val=Ryacas::yac(paste0("Expand(",tmp.val,")"),rettype="str")
          all.F3.eqs=c(all.F3.eqs,paste0(tmp.nom,tmp.val))
        }
      }
    }
  }
  #recoding to original branch length
  for(i in 1:x@n.edges){all.F3.eqs=gsub(tmp.edges.code[i],x@edges.names[i],all.F3.eqs)}
  #printing
  if(printfile){  
  cat(file=outfile,"\n\n#########################################\n############## F3 equations #############\n#########################################\n\n",append=TRUE)
  write(all.F3.eqs,file=outfile,append = TRUE)
  }
  
  #################
  #####F4
  ################    
  all.F4.eqs=c() #brute force since character (and fast enough (don't expect lots of pops!))
  for(i in 1:(x@n.leaves-1)){
    for(j in (i+1):x@n.leaves){
      tmp.popindex=(1:x@n.leaves)[-c(i,j)] 
      tmp.popindex=tmp.popindex[tmp.popindex>i] 
      tmp.popindex.lg=length(tmp.popindex)
      if(tmp.popindex.lg>1){
        for(k in 1:(tmp.popindex.lg-1)){
          for(l in (k+1):tmp.popindex.lg){  
            p=tmp.popindex[k] ; q=tmp.popindex[l] 
            tmp.nom=paste0("F4(",x@leaves[i],",",x@leaves[j],";",x@leaves[p],",",x@leaves[q],") = ")
            tmp.val=paste0("(",mat.f2[i,q],"+",mat.f2[j,p],"-(",mat.f2[i,p],"+",mat.f2[j,q],"))/2")
            tmp.val=Ryacas::yac(paste0("Expand(",tmp.val,")"),rettype="str")
            all.F4.eqs=c(all.F4.eqs,paste0(tmp.nom,tmp.val))
          }
        }
      }
    }
  }
  #recoding to original branch length
  for(i in 1:x@n.edges){all.F4.eqs=gsub(tmp.edges.code[i],x@edges.names[i],all.F4.eqs)}
  #printing
  if(printfile){  
  cat(file=outfile,"\n\n#########################################\n############## F4 equations #############\n#########################################\n\n",append=TRUE)
  write(all.F4.eqs,file=outfile,append = TRUE)
  }
  #################
  #####F2
  ################

  #recoding of mat.f2 to original branch length
 # for(i in 1:out@n.edges){mat.f2=gsub(tmp.edges.code[i],x@edges.names[i],mat.f2)}
  all.F2.eqs=c() #brute force since character (and fast enough (don't expect lots of pops!))
  for(i in 1:(x@n.leaves-1)){
    for(j in (i+1):x@n.leaves){
      tmp.nom=paste0("F2(",x@leaves[i],",",x@leaves[j],") = ")
      all.F2.eqs=c(all.F2.eqs,paste0(tmp.nom,mat.f2[i,j]))
    }
  }
  for(i in 1:x@n.edges){all.F2.eqs=gsub(tmp.edges.code[i],x@edges.names[i],all.F2.eqs)}
  #printing
  if(printfile){  
  cat(file=outfile,"\n\n#########################################\n############## F2 equations #############\n#########################################\n\n",append=TRUE)
  write(all.F2.eqs,file=outfile,append = TRUE)
  }
  
return(list(model.matrix=x.mat,omega=omega,F2.equations=all.F2.eqs,F3.equations=all.F3.eqs,F4.equations=all.F4.eqs))  
}
