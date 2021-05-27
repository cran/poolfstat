#' Estimate parameters of an admixture graph
#' @param graph.params An object of class graph.params containing graph information and relevant Fstats estimates (see the function generate.graph.params)
#' @param Q.lambda A scalar (usually small) to add to the diagonal elements of the error covariance matrix of fstats estimates (may improve numerical stability of its decomposition for large number of populations)
## #' @param LeastSquares If TRUE, a Least Square cost function is used (default=FALSE, see details)
#' @param eps.admix.prop A scalar defining admixture proportion domain (eps.admix.prop vary between eps.admix.prop and 1-eps.admix.prop)
#' @param edge.fact The multiplying factor of edges length in graph representation  
#' @param admix.fact The multiplying factor of admixture proportion in graph representation
#' @param compute.ci Derive 95\% Confidence Intervals for the parameters of the admixture graph (edge lengths and admixture rates)
#' @param drift.scaling If TRUE scale edge lengths in drift units (require estimates of leave heterozygosities)
#' @param outfileprefix The prefix of the dot file that will represent the graph (with extension ".dot"). If NULL, no graph file generated
#' @param verbose If TRUE extra information is printed on the terminal
#' @details Let \eqn{f} represent the n-length vector of basis target (i.e., observed)  F2 and F3 statistics and \eqn{g(e;a)=X(a)*e} the vector of their expected values given the vector of graph edges lengths \eqn{e} and the incidence matrix \eqn{X(a)} that depends on the structure of the graph and the admixture rates \eqn{a} (if there is no admixture in the graph, \eqn{X(a)} only contains 0 or 1). The function attempts to find the  \eqn{e} and \eqn{a} graph parameter values that minimize a cost (score of the model) defined as \eqn{S(e;a)=(f-g(e;a))'Q^{-1}(f-g(e;a))}. Assuming \eqn{f~N(g(e;a),Q)} (i.e., the observed f-statistics vector is multivariate normal distributed around an expected g vector specified by the admixture graph and a covariance structure empirically estimated), \eqn{S=-2log(L) - K} where \eqn{L} is the likelihood of the fitted graph and \eqn{K=n*log(2*pi)+log(|Q|)}. Also, for model comparison purpose, a standard \eqn{BIC} is then derived from \eqn{S} as \eqn{BIC= S + p*log(n) - K} (p being the number of graph parameters, i.e., edge lengths and admixture rates). 
#' As mentioned by Patterson et al. (2012), the score \eqn{S(e;a)} is quadratic in edge lengths \eqn{e} given \eqn{a}. The function uses the Lawson-Hanson non-negative linear least squares algorithm implemented in the nnls function (package nnls) to estimate \eqn{e} (subject to the constraint of positive edge lengths) by finding the vector \eqn{e} that minimize \eqn{S(e;a)=(f-X(a)*e)'Q^{-1}(f-X(a)*e)=||G*f-G*X(a)*e||} (where \eqn{G} results from the Cholesky decomposition of \eqn{Q^{-1}}, i.e., \eqn{Q^{-1}=G'G}). Note that the *Q.lambda* argument may be used to add a small constant (e.g., \eqn{1e-4}) to the diagonal elements of \eqn{Q} to avoid numerical problems (see Patterson et al., 2012). Yet *Q.lambda* is always disregarded when computing the final score \eqn{S} and \eqn{BIC}. Minimization of \eqn{S(e;a)} is thus reduced to the identification of the admixture rates (\eqn{a} vector) which is performed using the L-BFGS-B  method (i.e., Limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm with box constraints) implemented in the optim function (stats package). The *eps.admix.prop* argument allows specifying the lower and upper bound of the admixture rates to *eps.admix.prop* and *1-eps.admix.prop* respectively.
#' Scaling of the edges lengths in drift units (i.e., in units of \eqn{t/2N} where \eqn{t} is time in generations and \eqn{N} is the effective population size) is performed as described in Lipson et al. (MBE, 2013) by dividing the estimated edges lengths by half the estimated heterozygosity of their parental nodes (using the property \eqn{hp=hc+2e(C,P)} where \eqn{hp} and \eqn{hc} are the heterozygosities of a child C and its parent P node and \eqn{e(C,P)} is the estimated length of the branch relating C and P.
#' Finally, if compute.ci=TRUE, a (rough) \eqn{95\%} confidence intervals is computed using a bisection method (with a \eqn{1e-4} precision) for each parameters in turn (all others being set to their estimated value). Note that \eqn{95\%} CI are here defined as the set of values associated to a score \eqn{S} such that \eqn{Sopt<S<Sopt+3.84} (where \eqn{Sopt} is the optimized score), i.e., with a likelihood-ratio test statistic with respect to the fitted values \eqn{<3.84} (the \eqn{95\%} threshold of a one ddl Chi-square distribution).
#' @return An object of class fitted.graph (see help(fitted.graph) for details)
#' @seealso To generate a graph.params object, see \code{\link{generate.graph.params}}. The fitted graph may be plotted directly using plot that calls grViz() function and the resulting fitted fstats may be compared to the estimated ones with \code{\link{compare.fitted.fstats}}.
#' @export
fit.graph<-function(graph.params,Q.lambda=0,eps.admix.prop=1e-6,edge.fact=1000,admix.fact=100,compute.ci=F,drift.scaling=F,outfileprefix=NULL,verbose=TRUE){
  time1=proc.time()    
  if(!is.graph.params(graph.params)){stop("The input graph.params is not a valid graph.params object (see the function generate.graph.params)\n")}
  if(length(graph.params@f2.target)==0){stop("No target fstats found in the graph.params object (generate.graph.params should be run with fstats results)\n")}
  if(sum(is.na(Q.lambda))>0 | Q.lambda<0 | length(Q.lambda)!=1){
    cat("Warning: Q.lambda must be a positive real\nIt has been set to zero\n")
    Q.lambda=0
  }
  invQmat=graph.params@f.Qmat
  diag(invQmat)=diag(invQmat)+Q.lambda
  invQmat<-solve(invQmat)
  Gamma=chol(invQmat)
  
  npop=graph.params@n.leaves
  n.f2=nrow(graph.params@f2.target.pops) ;  n.f3=nrow(graph.params@f3.target.pops)  
  f.target=c(graph.params@f2.target,graph.params@f3.target)
  om.target=Gamma%*%f.target
  n.fstats=length(f.target)
  n.edges=ncol(graph.params@graph.matrix) #allows excluding admixture edges if any
  Pmat=graph.params@graph.matrix
  
  ######info to efficiently compute x.mat from Pmat (by using omega matrix element)
  ###### the coefficients (relative to the edge matrix) of the wXY elements of the omega matrix are obtained from Pmat
  ###### as follows
  ######  i)  for the npop diagonal elements Pmat*Pmat gives a npop by nedges matrix of edges weights
  ######  ii) for the npop*(npop-1)/2 off diagonal elements a npop*(npop-1)/2 x nedges matrix of edges weights (named coef.omega.offdiag) is obtained by taking the indexes of the each pairs of pop in Pmat (omega.path.pop.idx). The mat.idx.coef.omega matrix stores the index of the population in matrix form to allow retrieving the relevant indexes in the coef matrix without worrying of the order of the comparisons
  ######  From these two matrix we can obtain edges weights for the target F2 and F3 of interest using:
  ######   F2(A,B)   = wAA + wBB - 2*wAB
  ######   F3(x;A,B) = wxx + wAB - wxA - wxB
  
  omega.path.pop.idx=t(combn(npop,2))
  mat.idx.coef.omega=matrix(0,npop,npop) #permet de gere les ordres differents entre matrices pour recuperer la bonne ligne des coefs: ATTENTION ordre de remplissable des index dans la matrice!!
  mat.idx.coef.omega[lower.tri(mat.idx.coef.omega)]=1:nrow(omega.path.pop.idx) #npop*(npop-1)/2!
  mat.idx.coef.omega[upper.tri(mat.idx.coef.omega)]=t(mat.idx.coef.omega)[upper.tri(mat.idx.coef.omega)]
  #Pmat row indx of the different populations involved in the target F2 and F3 
  tmp.names=rownames(Pmat)
  idx.pop.ref=which(tmp.names==graph.params@popref)
  idx.F2.pB=sapply(graph.params@f2.target.pops[,2],function(x){which(tmp.names==x)})
  idx.F3.pA=sapply(graph.params@f3.target.pops[,2],function(x){which(tmp.names==x)}) 
  idx.F3.pB=sapply(graph.params@f3.target.pops[,3],function(x){which(tmp.names==x)})
  idx.off.diag.F2.wAB.term=(idx.pop.ref-1)*npop + idx.F2.pB
  idx.off.diag.F3.wxA.term=(idx.pop.ref-1)*npop + idx.F3.pA
  idx.off.diag.F3.wxB.term=(idx.pop.ref-1)*npop + idx.F3.pB
  idx.off.diag.F3.wAB.term=(idx.F3.pA-1)*npop + idx.F3.pB 
  
  ##internal function to derive xmat from an evaluate Pmat
  Pmat_to_xmat_<-function(PP){
    coef.omega.diag=PP*PP
    coef.omega.offdiag=PP[omega.path.pop.idx[,1],]*PP[omega.path.pop.idx[,2],]  
    tmp.f2.coef=coef.omega.diag[rep(idx.pop.ref,n.f2),] + coef.omega.diag[idx.F2.pB,] -
      2*coef.omega.offdiag[mat.idx.coef.omega[idx.off.diag.F2.wAB.term],]
    tmp.f3.coef=coef.omega.diag[rep(idx.pop.ref,n.f3),] + 
      coef.omega.offdiag[mat.idx.coef.omega[idx.off.diag.F3.wAB.term],] -
      coef.omega.offdiag[mat.idx.coef.omega[idx.off.diag.F3.wxA.term],] -
      coef.omega.offdiag[mat.idx.coef.omega[idx.off.diag.F3.wxB.term],]  
    XX=rbind(tmp.f2.coef,tmp.f3.coef)
    #on ne garde qu'on des deux root edges (par construction meme poids sur les deux) pour que X'X soit full rank!
    #en mettant le root combine restant en premiere colonne  
    return(cbind(XX[,graph.params@root.edges.idx[1]],XX[,-graph.params@root.edges.idx]))
  }
  
  ##################
  ####Evaluation par calcul de x.mat:   
  ##################  
  
  
  if(!graph.params@is.admgraph){#pas besoin d'optimisation!
    if(verbose){
      cat("Estimation started (direct algebraic solution)\n")
    }
    suppressWarnings(Pmat.eval<-matrix(as.numeric(Pmat),nrow = npop,ncol=n.edges,dimnames = dimnames(Pmat)))
    x.mat=Pmat_to_xmat_(Pmat.eval)
    x.mat.rank<-Matrix::rankMatrix(x.mat)
    if(x.mat.rank<ncol(x.mat)){is.singular=TRUE}else{is.singular=FALSE}
    if(!is.singular){
      Pi=Gamma%*%x.mat
      #    tmp.res=solve(t(Pi)%*%Pi)%*%t(Pi)%*%om.target
      tmp.res=nnls(Pi,om.target)$x
      edges.length=rep(0,n.edges)
      edges.length[-graph.params@root.edges.idx]=tmp.res[-1]
      edges.length[graph.params@root.edges.idx]=tmp.res[1]/2
      tmp=om.target-Pi%*%tmp.res
      tmp.score=as.numeric(t(tmp)%*%tmp)
      tmp.opt.res=list(edges.length=edges.length,admix.prop=NA,f.values=x.mat%*%tmp.res,score=tmp.score)  
    }
  }else{
    if(eps.admix.prop<0 | eps.admix.prop>0.1){stop("ERROR: eps.admix.prop must be a small positive value\n")}
    if(n.fstats< (n.edges+graph.params@n.adm.nodes)){
      is.singular=TRUE
    }else{
     admix.prop.inits=rep(0.5,graph.params@n.adm.nodes)
     Pmat.params<-unique(as.vector(Pmat))
     Pmat.vect<-match(Pmat,Pmat.params) #pour evaluer les 0 et les 1 qu'une seule fois 
     #test singularite: en mettant les coef a 0.5
     for(i in 1:graph.params@n.adm.nodes){assign(graph.params@adm.params.names[i],admix.prop.inits[i])}
     Pmat.params.eval=c()
     for(i in 1:length(Pmat.params)){Pmat.params.eval=c(Pmat.params.eval,eval(parse(text=Pmat.params[i])))}
     Pmat.eval=matrix(Pmat.params.eval[Pmat.vect],npop,n.edges) 
     x.mat=Pmat_to_xmat_(Pmat.eval)
     x.mat.rank<-Matrix::rankMatrix(x.mat)
     if(x.mat.rank<ncol(x.mat)){is.singular=TRUE}else{is.singular=FALSE}
    } 
    #########################  
    if(!is.singular){
      ff_optim_<-function(aa,return.score=TRUE){#si true retunr onlys score => pour optimisation, sinon tout pour le reste
        for(i in 1:graph.params@n.adm.nodes){assign(graph.params@adm.params.names[i],aa[i])}
        Pmat.params.eval=c()
        for(i in 1:length(Pmat.params)){Pmat.params.eval=c(Pmat.params.eval,eval(parse(text=Pmat.params[i])))}
        Pmat.eval=matrix(Pmat.params.eval[Pmat.vect],npop,n.edges) 
        x.mat=Pmat_to_xmat_(Pmat.eval)   
        Pi=Gamma%*%x.mat
        #      tmp.edges.length=solve(t(Pi)%*%Pi)%*%t(Pi)%*%om
        tmp.edges.length=nnls(Pi,om.target)$x
        tmp=om.target-Pi%*%tmp.edges.length
        tmp.score=as.numeric(t(tmp)%*%tmp)
        if(return.score){
          return(tmp.score)
        }else{
          edges.length=rep(0,n.edges)
          edges.length[-graph.params@root.edges.idx]=tmp.edges.length[-1]
          edges.length[graph.params@root.edges.idx]=tmp.edges.length[1]/2            
          return(list(edges.length=edges.length,admix.prop=aa,f.values=x.mat%*%tmp.edges.length,score=tmp.score))
        }
      }
      #  inits.par=c(-1*log((1-admix.prop.inits)/admix.prop.inits))
      if(verbose){
        cat("Starting estimation of admixture rates (LBFGS score optimisation)\n")
        # cat("\tInitial score=",ff_optim_(inits.par),"\n")
        cat("\tInitial score=",ff_optim_(admix.prop.inits),"\n")      
      }
      #    tmp.opt=lbfgs(ff_optim_, function(x) pracma::grad(ff_optim_, x),inits.par,invisible=1,max_iterations = 1000*length(inits.par)) 
      tmp.opt=optim(admix.prop.inits,ff_optim_, method = "L-BFGS-B",
                    lower = rep(eps.admix.prop,graph.params@n.adm.nodes),
                    upper = rep(1-eps.admix.prop,graph.params@n.adm.nodes),
                    control = list(maxit = 10000*graph.params@n.adm.nodes))
      tmp.opt.res=ff_optim_(tmp.opt$par,return.score = FALSE)
    }
  }
  
  if(is.singular){
    if(verbose){cat("The system is singular: the rank of the incidence matrix is lower than n.edges-1 (e.g., some branches are not identifiable) and/or the number of parameters is larger than the number of basis f-statistics\nIt is not possible to properly fit this admixture graph.\n")}
    out=NULL
  }else{
    #####################
    #####Calcul BIC + impression infos si verbose
    #####################
    
    nparams=graph.params@n.adm.nodes+graph.params@n.edges-1
    Kvalue=n.fstats*log(2*pi)+determinant(graph.params@f.Qmat,logarithm=TRUE)$modulus[1] #-2*cte de likelihood
    if(Q.lambda>0){#on repart de la bonne matrice pour caluler les vraisemblance
      invQmat=solve(graph.params@f.Qmat)
      final.score=f.target-tmp.opt.res$f.values
      final.score=as.numeric(t(final.score)%*%invQmat%*%final.score)
    }else{
      final.score=tmp.opt.res$score
    }
    bic=final.score+nparams*log(n.fstats) - Kvalue
    
    if(verbose){
      time.elapsed=(proc.time()-time1)[3]
      nminutes=floor(time.elapsed/60) ;  nseconds=round(time.elapsed-nminutes*60)
      cat("Estimation ended in",nminutes, "m ",nseconds,"s\n")
      cat("\tFinal Score:",final.score,"\n")
      cat("\tBIC:",bic,"\n")    
      if(graph.params@is.admgraph){
        if(tmp.opt$convergence!=0){
          cat("\nWARNING: possible convergence issues in optimisation\nSee out@optim.results$convergence and out@optim.results$message for details\n")}
      }    
    }    
    
    edges.length=tmp.opt.res$edges.length
    names(edges.length)=graph.params@edges.names
    admix_prop=tmp.opt.res$admix.prop
    fitted.outstats=cbind(f.target,tmp.opt.res$f.values,(tmp.opt.res$f.values-f.target)/sqrt(diag(graph.params@f.Qmat)))
    colnames(fitted.outstats)=c("Stat. value","Fitted Value","Z-score")
    rownames(fitted.outstats)=c(rownames(graph.params@f2.target.pops),rownames(graph.params@f3.target.pops))
    f2.mat=matrix(0,npop,npop,dimnames=list(graph.params@leaves,graph.params@leaves))
    tmp.f2=tmp.opt.res$f.values[(1:nrow(graph.params@f2.target.pops))]
    names(tmp.f2)=rownames(graph.params@f2.target.pops)
    tmp.f3=fitted.outstats[-(1:nrow(graph.params@f2.target.pops)),2]
    f2.mat[graph.params@f2.target.pops[1,1],graph.params@f2.target.pops[,2]]=
      f2.mat[graph.params@f2.target.pops[,2],graph.params@f2.target.pops[1,1]]=tmp.f2
    for(i in 1:length(tmp.f3)){
      tmp.p1=graph.params@f3.target.pops[i,2] ; tmp.p2=graph.params@f3.target.pops[i,3]
      f2.mat[tmp.p1,tmp.p2]=f2.mat[tmp.p2,tmp.p1]=tmp.f2[paste0(graph.params@popref,",",tmp.p1)] + tmp.f2[paste0(graph.params@popref,",",tmp.p2)] - 2*tmp.f3[i]
    }
    
    #################
    ####Computing CI (bisection method)
    #################
    if(compute.ci){
      tmp.edges.length=c(sum(edges.length[graph.params@root.edges.idx]),edges.length[-graph.params@root.edges.idx]) #on concatene les deux root edges
      like_eval_<-function(br.l=tmp.edges.length,aa=admix_prop){
        if(graph.params@is.admgraph){#sinon x.mat est deja evalue
          for(i in 1:graph.params@n.adm.nodes){assign(graph.params@adm.params.names[i],aa[i])}
          Pmat.params.eval=c()
          for(i in 1:length(Pmat.params)){Pmat.params.eval=c(Pmat.params.eval,eval(parse(text=Pmat.params[i])))}
          Pmat.eval=matrix(Pmat.params.eval[Pmat.vect],npop,n.edges) 
          x.mat=Pmat_to_xmat_(Pmat.eval)   
          tmp.outstats=x.mat%*%br.l
        }
        tmp=c(f.target-x.mat%*%br.l)
        return(as.numeric(t(tmp)%*%invQmat%*%tmp))
      }    
      #grid search
      #  i=1
      #   br.l.eval=seq(0.0017,0.0023,0.00001)
      #    tt=lapply(1:length(br.l.eval),nn<-function(x){tt=edges.length;tt[i]=br.l.eval[x];return(like_eval_(br.l=tt))})
      #    plot(br.l.eval,unlist(tt))
      #   abline(h=min(unlist(tt))+3.84)
      
      target.score=final.score+3.84 #95% ci
      eps.tol=1e-4
      
      tmp.br.l.ci=matrix(0,graph.params@n.edges-1,2) #avec root en premier
      for(i in 1:(graph.params@n.edges-1)){
        br.l.eval=tmp.edges.length
        #a droite
        tmp.binf=tmp.edges.length[i] ;  tmp.bsup=1.
        br.l.eval[i]=tmp.bsup ;   sup.score=like_eval_(br.l=br.l.eval)
        if(sup.score>target.score){
          while((tmp.bsup-tmp.binf)>eps.tol){
            tmp.cdt=mean(c(tmp.bsup,tmp.binf))
            br.l.eval[i]=tmp.cdt
            tmp.score.eval=like_eval_(br.l=br.l.eval)
            if(tmp.score.eval>target.score){tmp.bsup=tmp.cdt}else{tmp.binf=tmp.cdt}
          }
        }
        tmp.br.l.ci[i,2]=tmp.bsup
        #a gauche
        tmp.bsup=tmp.edges.length[i] ; tmp.binf=0.
        br.l.eval[i]=tmp.binf ;   inf.score=like_eval_(br.l=br.l.eval)
        if(inf.score>target.score){
          while((tmp.bsup-tmp.binf)>eps.tol){
            tmp.cdt=mean(c(tmp.bsup,tmp.binf))
            br.l.eval[i]=tmp.cdt
            tmp.score.eval=like_eval_(br.l=br.l.eval)
            if(tmp.score.eval>target.score){tmp.binf=tmp.cdt}else{tmp.bsup=tmp.cdt}
          }
        }
        tmp.br.l.ci[i,1]=tmp.binf
      }
      br.l.ci=matrix(0,graph.params@n.edges,2)
      rownames(br.l.ci)=c(graph.params@edges.names)
      colnames(br.l.ci)=c("95% Inf.","95% Sup.")
      br.l.ci[-graph.params@root.edges.idx,]=tmp.br.l.ci[-1,]
      br.l.ci[graph.params@root.edges.idx,]=matrix(rep(tmp.br.l.ci[1,]/2,2),nrow=2,byrow=T)
      #admixture params
      if(graph.params@is.admgraph){ 
        adm.prop.ci=matrix(0,graph.params@n.adm.nodes,2)
        rownames(adm.prop.ci)=graph.params@adm.params.names
        colnames(adm.prop.ci)=c("95% Inf.","95% Sup.")
        for(i in 1:graph.params@n.adm.nodes){
          adm.eval=admix_prop
          #a droite
          tmp.binf=admix_prop[i] ;  tmp.bsup=1.- eps.admix.prop
          adm.eval[i]=tmp.bsup ;   sup.score=like_eval_(aa=adm.eval)
          if(sup.score>target.score){
            while((tmp.bsup-tmp.binf)>eps.tol){
              tmp.cdt=mean(c(tmp.bsup,tmp.binf))
              adm.eval[i]=tmp.cdt
              tmp.score.eval=like_eval_(aa=adm.eval)
              if(tmp.score.eval>target.score){tmp.bsup=tmp.cdt}else{tmp.binf=tmp.cdt}
            }
          }
          adm.prop.ci[i,2]=tmp.bsup
          #a gauche
          tmp.bsup=admix_prop[i] ; tmp.binf=eps.admix.prop
          adm.eval[i]=tmp.binf ;   inf.score=like_eval_(aa=adm.eval)
          if(inf.score>target.score){
            while((tmp.bsup-tmp.binf)>eps.tol){
              tmp.cdt=mean(c(tmp.bsup,tmp.binf))
              adm.eval[i]=tmp.cdt
              tmp.score.eval=like_eval_(aa=adm.eval)
              if(tmp.score.eval>target.score){tmp.binf=tmp.cdt}else{tmp.bsup=tmp.cdt}
            }
          }
          adm.prop.ci[i,1]=tmp.binf
        }
      }
      
    }  
    
    #################
    ###Drift scaling
    #################
    
    if(drift.scaling){
      if(is.null(graph.params@Het)){
        cat("Warning: drift.scaling requires estimates of the leave heterozygosity (cf. computefstats and generate.graph.params functions\n")
        cat("No scaling will be done\n")
        drift.scaling=FALSE
      }else{
        if(compute.ci){
          br.l.ci=cbind(br.l.ci,br.l.ci)
          colnames(br.l.ci)=c("95% Inf.","95% Sup.","95% Inf. (drift scaled)","95% Sup. (drift scaled)")
        }
        tmp.edges=matrix(unlist(strsplit(names(edges.length),split="<->")),nrow=length(edges.length),byrow=T)
        tmp.nodes=unique(as.vector(tmp.edges))
        tmp.het=rep(NA,length(tmp.nodes))
        names(tmp.het)=tmp.nodes
        tmp.het[names(graph.params@Het)]=graph.params@Het
        for(i in names(graph.params@Het)){
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
          cat("Warning: Problem in estimating internal node heterozygosity\n")
          cat("No scaling will be done\n")
          drift.scaling=FALSE
        }else{
          edges.length.scaled=2*edges.length/tmp.het[tmp.edges[,1]]
          if(compute.ci){br.l.ci[,3:4]=2*br.l.ci[,3:4]/tmp.het[tmp.edges[,1]]}
          nodes.het=tmp.het
        }
      }
    }
    
    #################
    ####graph representation of the results
    #################
    
    ##header of the digraph file
    leaves.color="green";size="7.5,10" #for dot file
    outlines=c("digraph G {\n",paste0("size = \"",size,"\" ;\n"),"labelloc = \"t\" ;\n","")
    ##info leaves
    leaves.color=rep(leaves.color,graph.params@n.leaves)
    outlines=c(outlines,paste0("\"",graph.params@leaves,"\" [ label = \"",graph.params@leaves,"\",color=",leaves.color,",style=filled ] ;\n"),"")
    ##info edges
    tmp.nodes=matrix(unlist(strsplit(names(edges.length),split="<->")),ncol=2,byrow = T)
    if(drift.scaling){
      tmp.values=as.character(round(edges.length.scaled*edge.fact))
    }else{
      tmp.values=as.character(round(edges.length*edge.fact))
    }
    tmp.info=rep("style=plain",nrow(tmp.nodes))
    if(graph.params@is.admgraph){
      adm.graph.rows=nchar(graph.params@graph[,3])>0
      tmp.nodes=rbind(tmp.nodes,graph.params@graph[adm.graph.rows,2:1])
      tmp.info=c(tmp.info,rep("style=dotted",sum(adm.graph.rows)))
      for(i in 1:graph.params@n.adm.nodes){assign(graph.params@adm.params.names[i],admix_prop[i])}
      tmp.adm.eval=graph.params@graph[adm.graph.rows,3]
      for(i in 1:length(tmp.adm.eval)){
        dum.eval=eval(parse(text=tmp.adm.eval[i]))
        tmp.adm.eval[i]=round(admix.fact*dum.eval,digits=0)
        if(admix.fact==100){tmp.adm.eval[i]=paste0(tmp.adm.eval[i],"%")}
      }
      tmp.values=c(tmp.values,tmp.adm.eval)
    }
    tmp=paste0("\"",tmp.nodes[,1],"\"->\"",tmp.nodes[,2],"\" [ label=\"",tmp.values,"\",",tmp.info," ] ;\n")
    outlines=c(outlines,tmp,"}")
    if(!is.null(outfileprefix)){
      writeLines(outlines,con=paste0(outfileprefix,".dot"))
      cat("Graph input file in dot format written in",paste0(outfileprefix,".inputgraph.dot"),".\nIt can be visualized using the grViz function from the DiagrammeR package\nby running the command:",paste0("grViz(\"",outfileprefix,".inputgraph.dot\")"),"\n")
    }
    
    out<-new("fitted.graph")    
    out@edges.length=edges.length
    if(graph.params@is.admgraph){out@admix.prop=admix_prop}
    out@score=final.score ; out@fitted.outstats=fitted.outstats
    out@fitted.f2.mat=f2.mat
    out@dot.graph=outlines
    if(graph.params@is.admgraph){out@optim.results=tmp.opt}
    #   out@pval=pval
    out@bic=bic
    out@graph=graph.params@graph
    if(drift.scaling){
      out@edges.length.scaled=edges.length.scaled ;  out@nodes.het=nodes.het
    }
    if(compute.ci){
      if(graph.params@is.admgraph){out@admix.prop.ci=adm.prop.ci}
      out@edges.length.ci=br.l.ci
    }
  }  
  return(out)
}