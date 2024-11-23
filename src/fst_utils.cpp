#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

//' @title compute_snpQ1
//' @name compute_snpQ1
//' @rdname compute_snpQ1
//'
//' @description
//' Compute SNP-specific Q1 by averaging over all samples
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param weight Vector of length npop giving the weighting scheme (w=1 for allele count data and w=poolsize/(poolsize-1) for PoolSeq data)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute all the SNP-specific Q1 over all pop. samples (useful for Fst computation with method Identity). 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q1
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ1")]]
 Rcpp::NumericVector compute_snpQ1(Rcpp::NumericMatrix refcount,
                                   Rcpp::NumericMatrix totcount,
                                   Rcpp::NumericVector weight,
                                   Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int i,j,nsamp ;
   double curval ;
   
   Rcpp::NumericVector Q1val(nsnps) ; // (initialise 0)
   Progress p(nsnps, display_progress);
   
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     nsamp=0 ;
     for(j=0;j<npops;j++){
       if(totcount(i,j)>1){
         nsamp++ ;
         curval=(2.*refcount(i,j)*(totcount(i,j)-refcount(i,j)))/((totcount(i,j)-1.)*totcount(i,j)) ;
         curval*=weight(j) ; //it is actually H1=1-Q1
         Q1val(i)+=curval ;
       }
     }
     Q1val(i)=1.-Q1val(i)/nsamp ;
     p.increment();   
   }
   return Q1val ;
 }

//' @title compute_snpQ1rw
//' @name compute_snpQ1rw
//' @rdname compute_snpQ1rw
//'
//' @description
//' Compute SNP-specific Q1 over all samples using weighting averages of pop. Q1 (eq. A46 in Hivert et al., 2018)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param weight Vector of length npop giving the weighting scheme (w=1 for allele count data and w=poolsize/(poolsize-1) for PoolSeq data)
//' @param sampsize Vector of length npop giving the haploid sample size (not used for count data)
//' @param readcount Logical (if TRUE PoolSeq data assumed i.e. weights depending on haploid size, otherwise weights depend on total counts)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute all the SNP-specific Q1 over all pop. samples using weighting averages of pop. Q1 as in eq. A46 of Hivert et al., 2018 (useful for Fst computation with method Identity). 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q1
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ1rw")]]
Rcpp::NumericVector compute_snpQ1rw(Rcpp::NumericMatrix refcount,
                                    Rcpp::NumericMatrix totcount,
                                    Rcpp::NumericVector weight,
                                    Rcpp::NumericVector sampsize,
                                    Rcpp::LogicalVector readcount,
                                    Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   bool poolseq = is_true(any(readcount)) ;  
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int i,j ;
   double q1,numval,denval,dumw ;
   
   Rcpp::NumericVector Q1val(nsnps) ; // (initialise 0)
   Rcpp::NumericVector npairs(npops) ;
   Progress p(nsnps, display_progress);
   
   if(poolseq){npairs=sampsize*(sampsize-1.);}
   
   Q1val.fill(NA_REAL);
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     numval=0.;denval=0. ;
     for(j=0;j<npops;j++){
       if(totcount(i,j)>1){
         q1=(2.*refcount(i,j)*(totcount(i,j)-refcount(i,j)))/((totcount(i,j)-1.)*totcount(i,j)) ;//it is actually H1=1-Q1
         if(poolseq){
           q1=1.-q1*weight(j) ;
           dumw=npairs(j);
         }else{
           q1=1.-q1;
           dumw=(totcount(i,j)-1.)*totcount(i,j) ;
         }
         numval+= (q1*dumw);
         denval+= dumw;
       }
     }
     if(denval>0.){Q1val(i)=numval/denval ;}
     p.increment();   
   }
   return Q1val ;
 }

//' @title compute_snpQ2
//' @name compute_snpQ2
//' @rdname compute_snpQ2
//'
//' @description
//' Compute SNP-specific Q2 by averaging over all pairs of samples
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param pairs Matrix of npoppairsx2 giving the index for all the pairs of pops included in the computation
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute all the SNP-specific Q2 over all pop. pairs (useful for Fst computation with method Identity). 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q2
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ2")]]
 Rcpp::NumericVector compute_snpQ2(Rcpp::NumericMatrix refcount,
                                   Rcpp::NumericMatrix totcount,
                                   Rcpp::IntegerMatrix pairs,
                                   Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   int nsnps=refcount.nrow();
   int npairs=pairs.nrow();
   int i,j,p1,p2,nsamp ;
   double dumval ;
   Rcpp::NumericVector Q2val(nsnps) ; // (initialise 0)
   
   Progress p(nsnps, display_progress);
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     nsamp=0 ;
     for(j=0;j<npairs;j++){
       p1=pairs(j,0) ; p2=pairs(j,1) ; 
       dumval=totcount(i,p1)*totcount(i,p2) ;
       if(dumval>0){
         nsamp++ ;
         Q2val(i)+=(refcount(i,p1)*refcount(i,p2)+(totcount(i,p1)-refcount(i,p1))*(totcount(i,p2)-refcount(i,p2)))/dumval ;
       }
     }
     Q2val(i)/=nsamp ;
     p.increment();
   }
   return Q2val ;
 }

//' @title compute_snpQ2rw
//' @name compute_snpQ2rw
//' @rdname compute_snpQ2w
//'
//' @description
//' Compute SNP-specific Q2 by averaging over all pairs of samples using weighting averages of pairwise Q2 (eq. A47 in Hivert et al., 2018)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param pairs Matrix of npoppairsx2 giving the index for all the pairs of pops included in the computation
//' @param sampsize Vector of length npop giving the haploid sample size (not used for count data)
//' @param readcount Logical (if TRUE PoolSeq data assumed i.e. weights depending on haploid size, otherwise weights depend on total counts)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute SNP-specific Q2 by averaging over all pairs of samples using weighting averages of pairwise Q2 (eq. A47 in Hivert et al., 2018)
//' (useful for Fst computation with method Identity). 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q2
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ2rw")]]
 Rcpp::NumericVector compute_snpQ2rw(Rcpp::NumericMatrix refcount,
                                     Rcpp::NumericMatrix totcount,
                                     Rcpp::IntegerMatrix pairs,
                                     Rcpp::NumericVector sampsize,
                                     Rcpp::LogicalVector readcount,
                                     Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   bool poolseq = is_true(any(readcount)) ;  
   int nsnps=refcount.nrow();
   int npairs=pairs.nrow();
   int i,j,p1,p2 ;
   double q2,numval,denval,dumval ;
   Rcpp::NumericVector Q2val(nsnps) ; // (initialise 0)
   Rcpp::NumericVector ncomp(npairs) ; // (only useful for PoolSeq)
   

   if(poolseq){
     for(i=0;i<npairs;i++){
       p1=pairs(i,0) ; p2=pairs(i,1) ;
       ncomp[i]=sampsize(p1)*sampsize(p2);}
     }

   Q2val.fill(NA_REAL);
   Progress p(nsnps, display_progress);
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     numval=0.; denval=0. ;
     for(j=0;j<npairs;j++){
       p1=pairs(j,0) ; p2=pairs(j,1) ; 
       dumval=totcount(i,p1)*totcount(i,p2) ;
       if(dumval>0){
         q2=(refcount(i,p1)*refcount(i,p2)+(totcount(i,p1)-refcount(i,p1))*(totcount(i,p2)-refcount(i,p2))) ;
         if(poolseq){
           q2=q2/dumval ;
           numval+=(ncomp(j)*q2);
           denval+=ncomp(j) ;
         }else{
           numval+=q2 ;
           denval+=dumval;
         }
       }
     }
     if(denval>0.){Q2val(i)=numval/denval ;}
     p.increment();
   }
   return Q2val ;
 }


//' @title compute_snpHierFstAov
//' @name compute_snpHierFstAov
//' @rdname compute_snpHierFstAov
//'
//' @description
//' Compute SNP-specific MSI, MSP, MSG, nc, nc_p and nc_pp used to derived the Anova estimator of hier. Fst for allele count or read count data (Pool-Seq)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param hapsize Vector of length npop giving the haploid size of each pool (if one element <=0, counts are interpreted as count data)
//' @param popgrpidx Vector of length npop giving the index (coded from 0 to ngrp-1) of the group of origin
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute SNP-specific MSI, MSP, MSG, nc, nc_p and nc_pp used to derived the Anova estimator of hier. Fst for allele count or read count data (Pool-Seq)
//' 
//' @return Return a nsnpsx6 matrix with SNP-specific MSI, MSP, MSG, nc, nc_p and nc_pp
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpHierFstAov")]]
 Rcpp::NumericMatrix compute_snpHierFstAov(Rcpp::NumericMatrix refcount,
                                           Rcpp::NumericMatrix totcount,
                                           Rcpp::IntegerVector hapsize,
                                           Rcpp::IntegerVector popgrpidx,
                                           Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   bool poolseq = is_false(any(hapsize<=0)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int ngrps=Rcpp::max(popgrpidx)+1;//popgrpidx must be 0-indexed
   
   Rcpp::NumericMatrix outmat(nsnps,6); // (initialise 0); of MSI, MSP, MSG, nc, nc_p and nc_pp (in this column order)
   
   int i,j,g,npopsampled ;
   double dum,S1,S2,S3,SSI,SSP,SSG ;
   double Xooo_ref ; // sum counts over all pops / S1
   
   Rcpp::NumericVector S1g(ngrps),S2g(ngrps);
   Rcpp::NumericVector Xgoo_ref(ngrps); // sum ref counts of all pops within each group / S1g
   Rcpp::NumericVector Xgio_ref(npops); // sample freq in each pop (if c>0, pop ignored otherwise)
   
   //for poolseq data
   double D2,D2_p,D2_pp ;
   Rcpp::NumericVector hapsizefact(npops),Dterm(npops),D2g(ngrps);
   if(poolseq){
     for(j=0;j<npops;j++){
       hapsizefact(j)=(hapsize(j)-1.)/hapsize(j) ;
     }
   }
   
   //   Rcout << "The value of ngrps: " << ngrps << " " << S1g.length() << "\n";  
   Progress p(nsnps, display_progress);
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     npopsampled=0;
     S1g.fill(0.) ; S2g.fill(0.) ; Xooo_ref=0. ; Xgoo_ref.fill(0.);
     for(j=0;j<npops;j++){
       g=popgrpidx(j);
       if(totcount(i,j)>0){
         npopsampled++;
         S1g(g)+=totcount(i,j) ;
         S2g(g)+=(totcount(i,j)*totcount(i,j)) ;
         Xooo_ref+=refcount(i,j);
         Xgoo_ref(g)+=refcount(i,j);
         Xgio_ref(j)=(refcount(i,j)+0.)/totcount(i,j);
       }
     }
     if(is_true(any(S1g<0.5)) || (npopsampled==ngrps)){
       //no count for at least one group (or only one pop in each group)
       //=> the SNP is discarded and all output values are set to NA
       for(j=0;j<6;j++){outmat(i,j)=NA_REAL ;}
     }else{
       S1=sum(S1g) ; S2=sum(S2g) ; S3=0. ;
       Xooo_ref/=S1;
       for(g=0;g<ngrps;g++){
         S3+=S1g(g)*S1g(g); 
         Xgoo_ref(g)/=S1g(g);
       }
       //compute SSI, SSP and SSG
       SSI=0.;SSP=0.;
       for(j=0;j<npops;j++){
         if(totcount(i,j)>0){
           //SSI update
           dum=refcount(i,j)*(1.-Xgio_ref(j))*(1.-Xgio_ref(j)) ;
           dum+= (totcount(i,j)-refcount(i,j))*(0.-Xgio_ref(j))*(0.-Xgio_ref(j)) ;
           SSI+=dum ;
           //SSP update
           g=popgrpidx(j);
           dum=Xgio_ref(j)-Xgoo_ref(g);
           SSP+=totcount(i,j)*dum*dum;
         }
       }
       SSI=2.*SSI; SSP=2.*SSP ;
       //compute SSG
       SSG=0.;
       for(g=0;g<ngrps;g++){
         dum=(Xgoo_ref(g)-Xooo_ref);
         SSG+=S1g(g)*dum*dum;
       }
       SSG=2.*SSG ;
       
       if(poolseq){
         //compute poolsize correction
         D2=0.;D2_pp=0.;D2g.fill(0.);
         for(j=0;j<npops;j++){
           if(totcount(i,j)>0){
             g=popgrpidx(j);
             Dterm(j)=(totcount(i,j)+0.)/hapsize(j) + hapsizefact(j);
             D2+=Dterm(j);
             dum=Dterm(j)*totcount(i,j);
             D2g(g)+=dum;D2_pp+=dum;
           }
         }
         D2_pp/=S1;
         D2_p=0.;
         for(g=0;g<ngrps;g++){
           D2_p+=(D2g(g)/S1g(g)) ;
         }
         //compute MSG, MSP and MSG
         outmat(i,0)=SSI/(S1-D2);      //MSI
         outmat(i,1)=SSP/(D2-D2_p);    //MSP
         outmat(i,2)=SSG/(D2_p-D2_pp); //MSG     
         //compute nc, nc_p, nc_pp
         dum=0. ;
         for(g=0;g<ngrps;g++){dum+=S2g(g)/S1g(g);}
         outmat(i,3)=(S1-dum)/(D2-D2_p);            //nc
         outmat(i,4)=(S1 - S3/S1)/(D2_p - D2_pp);   //nc_p
         outmat(i,5)=(S2/S1 - dum)/(D2_p - D2_pp);  //nc_pp        
       }else{
         //compute MSG, MSP and MSG
         outmat(i,0)=SSI/(S1-npopsampled);     //MSI
         outmat(i,1)=SSP/(npopsampled-ngrps);  //MSP
         outmat(i,2)=SSG/(ngrps-1.);           //MSG
         //compute nc, nc_p, nc_pp
         dum=0. ;
         for(g=0;g<ngrps;g++){dum+=S2g(g)/S1g(g);}
         outmat(i,3)=(S1-dum)/(npopsampled-ngrps); //nc
         outmat(i,4)=(S1 - S3/S1)/(ngrps-1.);      //nc_p
         outmat(i,5)=(S2/S1 - dum)/(ngrps-1.);     //nc_pp
         //  Rcout << "SNP i: " << outmat(i,3) << " " << "\n"; 
       }
     }
     p.increment();
   }
   
   return outmat ;
 }


//' @title compute_snpFstAov
//' @name compute_snpFstAov
//' @rdname compute_snpFstAov
//'
//' @description
//' Compute SNP-specific MSG, MSP and nc used to derived the Anova estimator of Fst for allele count or read count data (Pool-Seq)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param hapsize Vector of length npop giving the haploid size of each pool (if one element <=0, counts are interpreted as count data)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute SNP-specific Q1 and Q2 based on Anova estimator of Fst for allele count or read count data (Pool-Seq).
//' For allele count data, the implemented estimator corresponds to that described in Weir, 1996 (eq. 5.2)  
//' For read (Pool-Seq) data, the implemented estimator corresponds to that described in Hivert et al., 2016  
//' 
//' @return Return a nsnpsx3 matrix with SNP-specific MSG, MSP and nc
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpFstAov")]]
 Rcpp::NumericMatrix compute_snpFstAov(Rcpp::NumericMatrix refcount,
                                       Rcpp::NumericMatrix totcount,
                                       Rcpp::IntegerVector hapsize,
                                       Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   bool poolseq = is_false(any(hapsize<=0)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   
   Rcpp::NumericMatrix outmat(nsnps,3) ; // (initialise 0); contains output values for all SNPs of Q1 and Q2 (in this column order)
   
   int i,j,npopsampled ;
   double dum,S1,S2,SSG,SSP,PA ;
   
   //for poolseq data
   double D2,D2star ;
   Rcpp::NumericVector hapsizefact(npops);
   if(poolseq){
     for(j=0;j<npops;j++){
       hapsizefact(j)=(hapsize(j)-1.)/hapsize(j) ;
     }
   }
   
   //   Rcout << "The value of ngrps: " << ngrps << " " << S1g.length() << "\n";  
   Progress p(nsnps, display_progress);
   for(i=0;i<nsnps;i++){
     if (Progress::check_abort() )
       return -1.0;
     npopsampled=0; S1=0.;S2=0.;PA=0.;
     //compute SumNi and npopsamples
     for(j=0;j<npops;j++){
       if(totcount(i,j)>0){
         npopsampled++;
         dum=(totcount(i,j)+0.) ;//casting counts
         S1+=dum ; S2+=(dum*dum) ;
         PA+=refcount(i,j);
       }
     }
     if(npopsampled<2){
       for(j=0;j<3;j++){outmat(i,j)=NA_REAL ;}
     }else{
       //compute SSG, SSP and PA (same for allele or read count)     
       PA/=S1;
       SSG=0.;SSP=0.;
       for(j=0;j<npops;j++){
         if(totcount(i,j)>0){
           dum=(refcount(i,j)+0.)/(totcount(i,j)+0.) ;//casting counts             
           SSG+=(refcount(i,j)*(1.- dum));
           dum-=PA;
           SSP+=(totcount(i,j)*dum*dum);
         }
       }
       SSG=2.*SSG ; SSP=2.*SSP ;
       //le facteur 2 dans le calcul de snp.Q1 et snp.Q2 vient du mode de calcul simplifie des MSP et MSG (en bi-allelique) avec les formules a la Weir
       //On peut retomber strictement sur les formules de Rousset 2007 (etendues dans Hivert et al. par Renaud pour le cas PoolSeq) qui matchent parfaitement l'approche poolseq ci-desous en faisant
       //S1 <- rowSums(x@total.count) ; S2 <- rowSums(x@total.count**2)
       //n_c=(S1-S2/S1)/(x@npops-1)
       //YY_alt=x@total.count - x@refallele.count
       //SSI = rowSums(x@refallele.count -  x@refallele.count^2 / x@total.count +  YY_alt - YY_alt^2 / x@total.count,na.rm = TRUE)
       //Ici on calcule les identites intra avec la formule:
       //  Q1=y - y^2/n + (n-y) - (n-y)^2/n
       //    =n - ((y+(n-y))^2 - 2y(n-y))/n
       //    =n - (n^2 - 2y(n-y))/n
       //    =2y(n-y)/n [utilise dans le calcul du MSG ci-dessus]
       //  SSP=rowSums(x@total.count * ((x@refallele.count / x@total.count) - (rowSums(x@refallele.count) / S1))^2 + x@total.count * ((YY_alt / x@total.count) - (rowSums(YY_alt) / S1))^2,na.rm = TRUE)
       //  MSI=SSI/(S1-x@npops) ; MSP=SSP/(x@npops-1)
       //  F_ST=(MSP-MSI)/(MSP+(n_c-1)*MSI)  #cf eq. 28A28 de Rousset (2007) en factorisant numerateur et denominateur par (S1-ns)*(ns-1) pour retrouver MSI et MSP
       //  F_ST_multi=mean(MSP-MSI,na.rm=T)/mean(MSP+(n_c-1)*MSI,na.rm=T)
       //  snp.Q1 = 1 - MSI  #cf eq 28A21 de Rousset 2007 en sommant sur tous les allele (i.e., SumPi_k=1)
       //  snp.Q2 = 1 - MSI - (MSP - MSI) / n_c     
       //      
       if(poolseq){
         //compute poolsize correction terms (D2 and D2star)
         D2=0. ; D2star=0. ;
         for(j=0;j<npops;j++){
           if(totcount(i,j)>0){
             dum=(totcount(i,j)+0.)/hapsize(j)+hapsizefact(j);
             D2+=dum ;
             D2star+=(totcount(i,j)*dum);
           }
         }
         D2star/=S1;
         outmat(i,0)=SSG/(S1 - D2); //MSG
         outmat(i,1)=SSP/(D2 - D2star);  //MSP 
         outmat(i,2)=(S1-(S2/S1))/(D2-D2star); //nc
       }else{
         outmat(i,0)=SSG/(S1-npopsampled) ; //MSG
         outmat(i,1)=SSP/(npopsampled-1.) ; //MSP
         outmat(i,2)=(S1-(S2/S1))/(npopsampled-1.); //nc
         //  Rcout << "SNP i: " << outmat(i,3) << " " << "\n"; 
       }
     }
     p.increment();
   }
   
   return outmat ;
 }

//' @title block_sum
//' @name block_sum
//' @rdname block_sum
//'
//' @description
//' Sugar to compute the sum of a stat per block
//'
//' @param stat vector of n stat values
//' @param snp_bj_id integer n-length vector with block index (from 0 to nblock-1) of the stat value 
//'
//' @details
//'  Sugar to compute the sum of a stat per block
//' 
//' @return Return a vector of length nblocks containing the per-block sums of the input stat
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".block_sum")]]
 Rcpp::NumericVector block_sum( Rcpp::NumericVector stat,
                                Rcpp::IntegerVector snp_bj_id){
   int nstats=stat.length();
   int nblocks=Rcpp::max(snp_bj_id) + 1 ;
   int i ;
   Rcpp::NumericVector bsum(nblocks) ; //initialise a 0
   for(i=0;i<nstats;i++){
     bsum(snp_bj_id(i))+=stat(i) ;
   }
   return bsum ;
 }

//' @title block_sum2
//' @name block_sum2
//' @rdname block_sum2
//'
//' @description
//' Sugar to compute the sum of a stat per block defined by a range of SNPs (allow treating overlapping blocks)
//'
//' @param stat vector of n stat values
//' @param snp_bj_id integer matrix of dim nblocks x 2 giving for each block the start and end stat value index 
//'
//' @details
//'  Sugar to compute the sum of a stat per block defined by a range of SNPs (allow treating overlapping blocks)
//' 
//' @return Return a vector of length nblocks containing the per-block sums of the input stat
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".block_sum2")]]
 Rcpp::NumericVector block_sum2( Rcpp::NumericVector stat,
                                 Rcpp::IntegerMatrix snp_bj_id){
   int nblocks=snp_bj_id.nrow() ;
   int i,j ;
   Rcpp::NumericVector bsum(nblocks) ; //initialise a 0
   for(i=0;i<nblocks;i++){
     //    Rcout << "Blk "<< i << " : " << snp_bj_id(i,0) << " " << snp_bj_id(i,1) << "\n"; 
     for(j=snp_bj_id(i,0);j<=snp_bj_id(i,1);j++){//important <= and not < : to cover all the range
       bsum(i)+=stat(j) ;
     }
   }
   return bsum ;
 }

