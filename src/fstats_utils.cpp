#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


//' @title compute_H1
//' @name compute_H1
//' @rdname compute_H1
//'
//' @description
//' Compute (uncorrected) 1-Q1 for each block-jackknife block (if any) and over all the SNPs (i.e., either within or outside blocks)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param nblocks Integer giving the number of block-jackknife blocs (may be 0 if no block-jackknife)
//' @param block_id Integer vector of length nsnps with the (0-indexed) id of the block to which each SNP belongs (-1 for SNPs outside blocks)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute all the (uncorrected) H1=1-Q1 for each block-jackknife block (if any) and overall SNPs (within or outside blocks). 
//' It is indeed more convenient to compute H1 (rather than Q1) to apply corrections afterwards within R function 
//' 
//' @return Return a matrix with npops rows and nblocks+1 column giving the mean H1 of each pop within each block and for all SNPs (last column)
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_H1")]]
Rcpp::NumericMatrix compute_H1(Rcpp::IntegerMatrix refcount,
                              Rcpp::IntegerMatrix totcount,
                              int nblocks,
                              Rcpp::IntegerVector block_id,
                              Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int i,j ;
   double dumval ;
   Rcpp::NumericMatrix H1val(npops,nblocks+1) ; // (initialise 0)
   Rcpp::NumericVector Nsnp(nblocks+1) ; // (initialise 0)   
   Progress p(npops, display_progress);
   
   for(i=0;i<npops;i++){
     if (Progress::check_abort() )
       return -1.0;
     for(j=0;j<=nblocks;j++){Nsnp(j)=0 ;}
     for(j=0;j<nsnps;j++){
        if(totcount(j,i)>1){
         Nsnp(nblocks)++ ;
         dumval=(2.*refcount(j,i)*(totcount(j,i)-refcount(j,i)))/((totcount(j,i)-1.)*totcount(j,i)) ;
         H1val(i,nblocks)+=dumval ;
         if(block_id(j)>=0){
           Nsnp(block_id(j))++ ;
           H1val(i,block_id(j))+=dumval ;
         }
       }
     }
     H1val(i,nblocks)/=Nsnp(nblocks) ;
     if(nblocks>0){
      for(j=0;j<nblocks;j++){
        H1val(i,j)/=Nsnp(j) ;
      }
     }
      p.increment();   
     }
   return H1val ;
 }

//' @title compute_Q2
//' @name compute_Q2
//' @rdname compute_Q2
//'
//' @description
//' Compute all Q2 for each block-jackknife block (if any) and overall SNPs (within or outside blocks)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param nblocks Integer giving the number of block-jackknife blocs (may be 0 if no block-jackknife)
//' @param block_id Integer vector of length nsnps with the (0-indexed) id of the block to which each SNP belongs (-1 for SNPs outside blocks)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute all Q2 for each block-jackknife block (if any) and overall SNPs (within or outside blocks). 
//' 
//' @return Return a matrix with npops*(npops-1)/2 and nblocks+1 column giving the mean Q2 of each pairwise pop comp. within each block and for all SNPs (last column)
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_Q2")]]
Rcpp::NumericMatrix compute_Q2(Rcpp::IntegerMatrix refcount,
                               Rcpp::IntegerMatrix totcount,
                               int nblocks,
                               Rcpp::IntegerVector block_id,
                               Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int i1,i2,j,cnt ;
   double dumval ;
   Rcpp::NumericMatrix Q2val(npops*(npops-1)/2,nblocks+1) ; // (initialise 0)
   Rcpp::NumericVector Nsnp(nblocks+1) ; // (initialise 0)   
   Progress p(npops*(npops-1)/2, display_progress);

   cnt=0 ;   
   for(i1=0;i1<(npops-1);i1++){
     for(i2=i1+1;i2<npops;i2++){
       if (Progress::check_abort() )
         return -1.0;
       for(j=0;j<=nblocks;j++){Nsnp(j)=0 ;}
       for(j=0;j<nsnps;j++){
       dumval=1.*totcount(j,i1)*totcount(j,i2) ;
       if(dumval>0){
         Nsnp(nblocks)++ ;
         dumval=(refcount(j,i1)*refcount(j,i2)+(totcount(j,i1)-refcount(j,i1))*(totcount(j,i2)-refcount(j,i2)))/dumval ;
         Q2val(cnt,nblocks)+=dumval ;
         if(block_id(j)>=0){
           Nsnp(block_id(j))++ ;
           Q2val(cnt,block_id(j))+=dumval ;
         }
       }
     }
      Q2val(cnt,nblocks)/=Nsnp(nblocks) ;
     if(nblocks>0){
       for(j=0;j<nblocks;j++){
         Q2val(cnt,j)/=Nsnp(j) ;
       }
     }
     cnt++ ;
     p.increment();
   }
    }
   return Q2val ;
 }


//' @title poppair_idx
//' @name poppair_idx
//' @rdname poppair_idx
//'
//' @description
//' Compute the index of the pairwise comparison from the idx of each pop
//'
//' @param idx_pop1 Integer giving the (0-indexed) index of the first pop 
//' @param idx_pop2 Integer giving the (0-indexed) index of the second pop 
//' @param nidx Integer giving the total number of indexes (i.e., number of pops)
//'
//' @details
//' If idx_pop2 < idx_pop1, indexes are reversed
//' 
//' @return Return the (0-indexed) index for the row associated to the pairwise comparison in the ordered flat list of all (npop*(npop-1))/2 pairwise stats
//' 
//' @examples
//' #
int poppair_idx(int idx_pop1,
                 int idx_pop2,
                 int nidx){
   int idx;
   if(idx_pop1<idx_pop2){
     idx=idx_pop1*nidx + (idx_pop2+1) - (idx_pop1+1)*(idx_pop1+2)/2 - 1 ;
    // ((tmp.cb[1,]-1)*npops + tmp.cb[2,])-(tmp.cb[1,]*(tmp.cb[1,]+1)/2)
   }else{
     idx=idx_pop2*nidx + (idx_pop1+1) - (idx_pop2+1)*(idx_pop2+2)/2 - 1 ; 
   }
   return idx;
 }
 
//' @title compute_F3fromF2
//' @name compute_F3fromF2
//' @rdname compute_F3fromF2
//'
//' @description
//' Compute all F3 from overall F2 values
//'
//' @param F2val Numeric vector of length nF2=(npop*(npop-1))/2 with all pairwise F2 estimates
//' @param Hval Numeric vector of length npop with all within pop heterozygosity estimates
//' @param npops Integer giving the number of populations
//'
//' @details
//' Compute F3 and F3star estimates from F2 (and heterozygosities)
//' 
//' @return Return a matrix of length nF3=npops*(npops-1)*(npops-2)/2 rows and 2 columns corresponding to the F3 and F3star estimates
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_F3fromF2")]]
 Rcpp::NumericMatrix compute_F3fromF2( Rcpp::NumericVector F2val,
                                       Rcpp::NumericVector Hval,
                                       int npops){
   //                                            Rcpp::IntegerMatrix f2idx
   int nf3=npops*(npops-1)*(npops-2)/2 ; //npop*choose(npop-1,2)
   int cnt=0, i, j , k ;
   Rcpp::NumericMatrix f3val(nf3,2) ; //F3 and F3star
   
   for(i=0;i<npops;i++){
     for(j=0;j<npops-1;j++){
       for(k=j+1;k<npops;k++){
         if(j!=i && k!=i){
           f3val(cnt,0)=(F2val(poppair_idx(i,j,npops)) + F2val(poppair_idx(i,k,npops)) - F2val(poppair_idx(j,k,npops)))/2. ;            
           f3val(cnt,1)=f3val(cnt,0)/Hval(i) ;
           cnt++ ;  
           }
         }
       }
     }
   return f3val ;
 } 
 
 
//' @title compute_F3fromF2samples
//' @name compute_F3fromF2samples
//' @rdname compute_F3fromF2samples
//'
//' @description
//' Compute all F3 from F2 values obtained from each block-jackknife bloc
//'
//' @param blockF2 Numeric Matrix with nF2=(npop*(npop-1))/2 rows and nblocks columns matrix containing pairwise-pop F2 estimates for each block-jackknife sample (l.o.o.)
//' @param blockHet Numeric Matrix with npop rows and nblocks columns containing all within pop heterozygosity estimates for each block-jackknife sample (l.o.o.)
//' @param npops Integer giving the number of populations
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute F3 and F3star estimates and their s.e. based on block-jackknife estimates of all F2 (and heterozygosities)
//' 
//' @return Return a matrix with nF3=npops*(npops-1)*(npops-2)/2 rows and four columns corresponding to the mean and the s.e. of F3 and the mean and s.e. of F3star
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_F3fromF2samples")]]
Rcpp::NumericMatrix compute_F3fromF2samples( Rcpp::NumericMatrix blockF2,
                                             Rcpp::NumericMatrix blockHet,
                                             int npops,   
                                             Rcpp::LogicalVector verbose){
//                                            Rcpp::IntegerMatrix f2idx
  bool display_progress = is_true(any(verbose)) ;
//  int npops=f2idx.nrow();
  int nblocs=blockF2.ncol() ;
  int nf3=npops*(npops-1)*(npops-2)/2 ; //npop*choose(npop-1,2)
  int cnt=0, i, j , k, b ;
  double sdcorr=sqrt(nblocs-1.) ; //Bessel correction x (nblocs-1)/sqrt(nblocs) = sqrt(n/(n-1))*(n-1)/sqrt(n);
  double dumval ;
  Rcpp::NumericMatrix f3val(nf3,4) ; //mean and s.e. (initialise 0)
  Progress p(nf3, display_progress);
  
  for(i=0;i<npops;i++){
    for(j=0;j<npops-1;j++){
      for(k=j+1;k<npops;k++){
        if(j!=i && k!=i){
          if (Progress::check_abort() )
            return -1.0;
          for(b=0;b<nblocs;b++){
 //           dumval=(blockF2(f2idx(i,j),b) + blockF2(f2idx(i,k),b) - blockF2(f2idx(j,k),b))/2. ;
            dumval=(blockF2(poppair_idx(i,j,npops),b) + blockF2(poppair_idx(i,k,npops),b) - blockF2(poppair_idx(j,k,npops),b))/2. ;            
            f3val(cnt,0)+=dumval ; f3val(cnt,1)+=dumval*dumval ;
            dumval/=blockHet(i,b) ;
            f3val(cnt,2)+=dumval ; f3val(cnt,3)+=dumval*dumval ;            
           }
          for(b=0;b<4;b++){
            f3val(cnt,b)/=nblocs ;
          }
          f3val(cnt,1)=sqrt((f3val(cnt,1) - (f3val(cnt,0)*f3val(cnt,0))))*sdcorr ;
          f3val(cnt,3)=sqrt((f3val(cnt,3) - (f3val(cnt,2)*f3val(cnt,2))))*sdcorr ;          
          cnt++ ;  
          p.increment(); 
  }
  }
      }
    }
  return f3val ;
}


//' @title generateF3names
//' @name generateF3names
//' @rdname generateF3names
//'
//' @description
//' Generate all names for F3 stats (same order as computation)
//'
//' @param popnames String vector with the names of all the pops
//'
//' @details
//' Generate all the npops*(npops-1)*(npops-2)/2 names for F3 stats (same order as computation)
//' 
//' @return Return a string matrix with 4 columns including the complete F3 configuration names (of the form Px;P1,P2), and the names of each pop involved in the configuration
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".generateF3names")]]
 Rcpp::StringMatrix generateF3names( Rcpp::StringVector popnames){
   int npops=popnames.size();
   int nf3=npops*(npops-1)*(npops-2)/2 ; //npop*choose(npop-1,2)
   int cnt=0, i, j , k;
   Rcpp::StringMatrix f3names(nf3,4) ;
   
   for(i=0;i<npops;i++){
     for(j=0;j<npops-1;j++){
       for(k=j+1;k<npops;k++){
         if(j!=i && k!=i){
           f3names(cnt,0)=popnames(i) ;f3names(cnt,0)+=";";
           f3names(cnt,0)+=popnames(j); f3names(cnt,0)+=","; f3names(cnt,0)+=popnames(k) ;
           f3names(cnt,1)=popnames(i) ; f3names(cnt,2)=popnames(j) ;f3names(cnt,3)=popnames(k) ; 
           cnt++ ;  
         }
       }
     }
   }
   return f3names ;
 }

//' @title compute_F4fromF2
//' @name compute_F4fromF2
//' @rdname compute_F4fromF2
//'
//' @description
//' Compute all F4 from overall F2 and Q2 values
//'
//' @param F2val Numeric vector of length nF2=(npop*(npop-1))/2 with all pairwise F2 estimates
//' @param npops Integer giving the number of populations
//'
//' @details
//' Compute F4 from F2 (and heterozygosities)
//' 
//' @return Return a vector of length nF4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 rows corresponding to all the F4 estimates for all possible configurations
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_F4fromF2")]]
 Rcpp::NumericVector compute_F4fromF2( Rcpp::NumericVector F2val,
                                        int npops){
   int nf4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 ; //choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
   int cnt=0, i, j , k, l ;
   Rcpp::NumericVector f4val(nf4) ; //F4 and D
   
   for(i=0;i<npops-1;i++){
     for(j=i+1;j<npops;j++){
       for(k=i+1;k<npops-1;k++){
         for(l=k+1;l<npops;l++){
           if(k!=j && l!=j){
             f4val(cnt)=  F2val(poppair_idx(i,l,npops)) + F2val(poppair_idx(j,k,npops)) ;
             f4val(cnt)-= F2val(poppair_idx(i,k,npops)) + F2val(poppair_idx(j,l,npops)) ;
             f4val(cnt)/=2. ;
             cnt++ ;  
           }
         }
       }}}
   return f4val ;
 }
   
//' @title compute_F4fromF2samples
//' @name compute_F4fromF2samples
//' @rdname compute_F4fromF2samples
//'
//' @description
//' Compute all F4 from F2 values obtained from each block-jackknife bloc
//'
//' @param blockF2 Numeric Matrix with nF2=(npop*(npop-1))/2 rows and nblocks columns matrix containing pairwise-pop F2 estimates for each block-jackknife sample (l.o.o.)
//' @param npops Integer giving the number of populations
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute F4 estimates and their s.e. based on block-jackknife estimates of all F2 (and heterozygosities)
//' 
//' @return Return a matrix with nF4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 rows and two columns corresponding to the mean and the s.e. of F4 estimates for all possible configurations
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_F4fromF2samples")]]
 Rcpp::NumericMatrix compute_F4fromF2samples( Rcpp::NumericMatrix blockF2,
                                              int npops,
                                              Rcpp::LogicalVector verbose){
  // Rcpp::IntegerMatrix f2idx,
   bool display_progress = is_true(any(verbose)) ;
   // int npops=f2idx.nrow();
   int nblocs=blockF2.ncol() ;
   int nf4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 ; //choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
   int cnt=0, i, j , k, l, b ;
   double sdcorr=sqrt(nblocs-1.) ; //Bessel correction x (nblocs-1)/sqrt(nblocs) = sqrt(n/(n-1))*(n-1)/sqrt(n);
   double dumval ;
   Rcpp::NumericMatrix f4val(nf4,2) ; //mean and s.e. (initialise 0)
   Progress p(nf4, display_progress);
   
   for(i=0;i<npops-1;i++){
     for(j=i+1;j<npops;j++){
       for(k=i+1;k<npops-1;k++){
         for(l=k+1;l<npops;l++){
           if(k!=j && l!=j){
             if (Progress::check_abort() )
               return -1.0;
             for(b=0;b<nblocs;b++){
//               dumval=(blockF2(f2idx(i,l),b) + blockF2(f2idx(j,k),b) - blockF2(f2idx(i,k),b) - blockF2(f2idx(j,l),b))/2. ;
               dumval=  blockF2(poppair_idx(i,l,npops),b) + blockF2(poppair_idx(j,k,npops),b) ;
               dumval-= blockF2(poppair_idx(i,k,npops),b) + blockF2(poppair_idx(j,l,npops),b) ;
               dumval/=2. ;
               f4val(cnt,0)+=dumval ; f4val(cnt,1)+=dumval*dumval ;
           }
             f4val(cnt,0)/=nblocs ; f4val(cnt,1)/=nblocs ;
             f4val(cnt,1)=sqrt((f4val(cnt,1) - (f4val(cnt,0)*f4val(cnt,0))))*sdcorr ;
             cnt++ ;  
             p.increment(); 
         }
       }
         }}}
   return f4val ;
 }

//' @title compute_F4DfromF2samples
//' @name compute_F4DfromF2samples
//' @rdname compute_F4DfromF2samples
//'
//' @description
//' Compute all F4 and Dstat from F2 values obtained from each block-jackknife bloc
//'
//' @param blockF2 Numeric Matrix with nF2=(npop*(npop-1))/2 rows and nblocks columns matrix containing pairwise-pop F2 estimates for each block-jackknife sample (l.o.o.)
//' @param blockDenom Numeric Matrix with nF4=(npops*(npops-1)/2)*((npops-2)*(npops-3)/2)/2 rows and nblocks containing the estimates of the denominator of Dstat (see compute_blockDdenom) for each block-jackknife sample (l.o.o.) 
//' @param npops Integer giving the number of populations
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute F4 and D estimates and their s.e. based on block-jackknife estimates of all F2 (and heterozygosities)
//' 
//' @return Return a matrix with nF4=(npops*(npops-1)/2)*((npops-2)*(npops-3)/2)/2 rows and four columns corresponding to the mean and the s.e. of F4 and the mean and s.e. of Dstat
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_F4DfromF2samples")]]
Rcpp::NumericMatrix compute_F4DfromF2samples(Rcpp::NumericMatrix blockF2,
                                             Rcpp::NumericMatrix blockDenom,
                                             int npops,
                                             Rcpp::LogicalVector verbose){
   // Rcpp::IntegerMatrix f2idx,
   bool display_progress = is_true(any(verbose)) ;
   // int npops=f2idx.nrow();
   int nblocs=blockF2.ncol() ;
   int nf4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 ; //choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
   int cnt=0, i, j , k, l, b ;
   double sdcorr=sqrt(nblocs-1.) ; //Bessel correction x (nblocs-1)/sqrt(nblocs) = sqrt(n/(n-1))*(n-1)/sqrt(n);
   double dumval ;
   Rcpp::NumericMatrix f4val(nf4,4) ; //mean and s.e. of F4 and D (initialise 0)
   Progress p(nf4, display_progress);
   
   for(i=0;i<npops-1;i++){
     for(j=i+1;j<npops;j++){
       for(k=i+1;k<npops-1;k++){
         for(l=k+1;l<npops;l++){
           if(k!=j && l!=j){
             if (Progress::check_abort() )
               return -1.0;
             for(b=0;b<nblocs;b++){
               //               dumval=(blockF2(f2idx(i,l),b) + blockF2(f2idx(j,k),b) - blockF2(f2idx(i,k),b) - blockF2(f2idx(j,l),b))/2. ;
               dumval=  blockF2(poppair_idx(i,l,npops),b) + blockF2(poppair_idx(j,k,npops),b) ;
               dumval-= blockF2(poppair_idx(i,k,npops),b) + blockF2(poppair_idx(j,l,npops),b) ;
               dumval/=2. ;
               f4val(cnt,0)+=dumval ; f4val(cnt,1)+=dumval*dumval ;
               dumval/=blockDenom(cnt,b) ;
               f4val(cnt,2)+=dumval ; f4val(cnt,3)+=dumval*dumval ;               
             }
             for(b=0;b<4;b++){
               f4val(cnt,b)/=nblocs ;
             }
             f4val(cnt,1)=sqrt((f4val(cnt,1) - (f4val(cnt,0)*f4val(cnt,0))))*sdcorr ;
             f4val(cnt,3)=sqrt((f4val(cnt,3) - (f4val(cnt,2)*f4val(cnt,2))))*sdcorr ;  
             cnt++;
             p.increment(); 
           }
         }
       }}}
   return f4val ;
 }

//' @title compute_blockDdenom
//' @name compute_blockDdenom
//' @rdname compute_blockDdenom
//'
//' @description
//' Compute the denominator of the Dstat for all quadruplet configuration and each block-jackknife block (if any) and overall SNPs (within or outside blocks)
//'
//' @param refcount Matrix of nsnpxnpop with counts (genotype or reads) for the reference allele
//' @param totcount Matrix of nsnpxnpop with total counts or read coverages
//' @param nblocks Integer giving the number of block-jackknife blocs (may be 0 if no block-jackknife)
//' @param block_id Integer vector of length nsnps with the (0-indexed) id of the block to which each SNP belongs (-1 for SNPs outside blocks)
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute the denominator of the Dstat for all quadruplet configuration and each block-jackknife block (if any) and overall SNPs (within or outside blocks)
//' 
//' @return Return a matrix with nf4=(npops*(npops-1)/2)*((npops-2)*(npops-3)/2)/2 rows and nblocks+1 columns giving the mean Dstat-denominator (1-Q2ab)(1-Q2cd)
//'  for all quadruplet configuration and within each block-jackknife sample and over all SNPs (last column)
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_blockDdenom")]]
Rcpp::NumericMatrix compute_blockDdenom(Rcpp::IntegerMatrix refcount,
                                        Rcpp::IntegerMatrix totcount,
                                        int nblocks,
                                       Rcpp::IntegerVector block_id,
                                       Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   int nsnps=refcount.nrow();
   int npops=refcount.ncol();  
   int nf4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 ; //choose(npop,2)*choose(npop-2,2)/2 : on div par pour virer permutations (P1,P2;P3,P4)=(P3,P4;P1,P2)
   int cnt=0, i, j , k, l, b,s ;
   double dumval,q2AB,q2CD ;
   Rcpp::NumericMatrix Ddenom(nf4,nblocks+1) ; // (initialise 0)
   Rcpp::NumericVector Nsnp(nblocks+1) ; // (initialise 0)   
   Progress p(nf4, display_progress);
   
   cnt=0 ;  
   for(i=0;i<npops-1;i++){
     for(j=i+1;j<npops;j++){
       for(k=i+1;k<npops-1;k++){
         for(l=k+1;l<npops;l++){
           if(k!=j && l!=j){
             if (Progress::check_abort() )
               return -1.0;
             for(b=0;b<=nblocks;b++){Nsnp(b)=0 ;}
             for(s=0;s<nsnps;s++){
               q2AB=1.*totcount(s,i)*totcount(s,j) ; 
               q2CD=1.*totcount(s,k)*totcount(s,l) ;                
               if(q2AB>0 && q2CD>0){
                 Nsnp(nblocks)++ ;
                 q2AB=(refcount(s,i)*refcount(s,j)+(totcount(s,i)-refcount(s,i))*(totcount(s,j)-refcount(s,j)))/q2AB ;
                 q2CD=(refcount(s,k)*refcount(s,l)+(totcount(s,k)-refcount(s,k))*(totcount(s,l)-refcount(s,l)))/q2CD ;                 
                 dumval=(1.-q2AB)*(1.-q2CD) ;
                 Ddenom(cnt,nblocks)+=dumval ;
                 if(block_id(s)>=0){
                   Nsnp(block_id(s))++ ;
                   Ddenom(cnt,block_id(s))+=dumval ;
                 }
               }
             }
             Ddenom(cnt,nblocks)/=Nsnp(nblocks) ;
             if(nblocks>0){
               for(b=0;b<nblocks;b++){
                 Ddenom(cnt,b)/=Nsnp(b) ;
               }
             }
             cnt++ ;
             p.increment();
     }
           }}}}

   return Ddenom ;
 }

//' @title generateF4names
//' @name generateF4names
//' @rdname generateF4names
//'
//' @description
//' Generate all names for F4 stats (same order as computation)
//'
//' @param popnames String vector with the names of all the pops
//'
//' @details
//' Generate all the nf4=(npops*(npops-1)/2)*((npops-2)*(npops-3)/2)/2 names for F4 stats (same order as computation)
//' 
//' @return Return a string matrix with 5 columns including the complete F4 configuration names (of the form P1,P2;P3,P4), and the names of each pop involved in the configuration
//' 
//' #
//' @export
// [[Rcpp::export(name=".generateF4names")]]
Rcpp::StringMatrix generateF4names( Rcpp::StringVector popnames){
   int npops=popnames.size();
   int nf4=(npops*(npops-1)/2) * ((npops-2)*(npops-3)/2) / 2 ; //npop*choose(npop-1,2)
   int cnt=0, i, j , k,l;
   Rcpp::StringMatrix f4names(nf4,5) ;
   
   for(i=0;i<npops-1;i++){
     for(j=i+1;j<npops;j++){
       for(k=i+1;k<npops-1;k++){
         for(l=k+1;l<npops;l++){
           if(k!=j && l!=j){
           f4names(cnt,0)=popnames(i) ;f4names(cnt,0)+=","; f4names(cnt,0)+=popnames(j) ;f4names(cnt,0)+=";";
           f4names(cnt,0)+=popnames(k) ;f4names(cnt,0)+=","; f4names(cnt,0)+=popnames(l) ;
           f4names(cnt,1)=popnames(i) ; f4names(cnt,2)=popnames(j) ;f4names(cnt,3)=popnames(k) ;f4names(cnt,4)=popnames(l) ; 
           cnt++ ;  
         }
       }
     }
   }
    }
   return f4names ;
 }

//' @title bjack_cov
//' @name bjack_cov
//' @rdname bjack_cov
//'
//' @description
//' Compute the block-jackknife covariance between two stats
//'
//' @param stat1 Vector of block-jackknife values for the first stat
//' @param stat2 Vector of block-jackknife values for the second stat
//'
//' @details
//'  Compute the block-jackknife covariance between two stats with correction
//' 
//' @return Covariance values
//' 
//' @examples
//' #
double bjack_cov(Rcpp::NumericVector stat1,
                 Rcpp::NumericVector stat2){
   int nb=stat1.size();
   int i ;
   double sum_xy=0.,sum_x=0.,sum_y=0.;
   for(i=0;i<nb;i++){
     sum_xy+=stat1(i)*stat2(i) ;
     sum_x+=stat1(i) ; sum_y+=stat2(i) ; 
   }
   sum_xy/=nb ; sum_x/=nb ; sum_y/=nb ;
   
   return (sum_xy-sum_x*sum_y)*(nb-1.);// hat_cov=(sum_xy-sum_x*sum_y)*(nb/(nb-1)) and bj_cov=hat_cov*(nb-1)*(nb-1)/nb
 }


//' @title compute_QmatfromF2samples
//' @name compute_QmatfromF2samples
//' @rdname compute_QmatfromF2samples
//'
//' @description
//' Compute the Qmat matrix (error covariance between all F2 and F3 measures) from F2 block-jackknife estimates
//'
//' @param blockF2 Numeric Matrix with nF2=(npop*(npop-1))/2 rows and nblocks columns matrix containing pairwise-pop F2 estimates for each block-jackknife sample (l.o.o.)
//' @param npops Integer giving the number of populations
//' @param verbose Logical (if TRUE progression bar is printed on the terminal)
//'
//' @details
//' Compute the error covariance matrix Qmat (between all F2 and F3 measures) from F2 block-jackknife estimates (by recomuting all F3 for all blocks)
//' 
//' @return Return the (nF2+nF3)*(nF2+nF3) error covariance (symmetric) matrix
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_QmatfromF2samples")]]
 Rcpp::NumericMatrix compute_QmatfromF2samples( Rcpp::NumericMatrix blockF2,
                                              int npops,   
                                              Rcpp::LogicalVector verbose){
   bool display_progress = is_true(any(verbose)) ;
   //  int npops=f2idx.nrow();
   int nblocs=blockF2.ncol() ;
   int nf2=blockF2.nrow() ;
   int nf3=npops*(npops-1)*(npops-2)/2 ; //npop*choose(npop-1,2)
   int nstats=nf2+nf3 ;
   int cnt=0, i, j , k, b ;
   Rcpp::NumericVector val1(nblocs),val2(nblocs) ;   
   Rcpp::NumericMatrix blockF3(nf3,nblocs) ; //F3 values
   Rcpp::NumericMatrix Qmat(nstats,nstats) ; //covariance
   Progress p(nf3+nstats*(nstats-1)/2, display_progress);
 // Computing block F3  
   for(i=0;i<npops;i++){
     for(j=0;j<npops-1;j++){
       for(k=j+1;k<npops;k++){
         if(j!=i && k!=i){
           if (Progress::check_abort() )
             return -1.0;
           for(b=0;b<nblocs;b++){
             //           dumval=(blockF2(f2idx(i,j),b) + blockF2(f2idx(i,k),b) - blockF2(f2idx(j,k),b))/2. ;
             blockF3(cnt,b)=(blockF2(poppair_idx(i,j,npops),b) + blockF2(poppair_idx(i,k,npops),b) - blockF2(poppair_idx(j,k,npops),b))/2. ;            
           }
           cnt++ ;  
           p.increment(); 
         }
       }
     }
   }
 // Computing covariance
 for(i=0;i<nstats-1;i++){
   if (Progress::check_abort() )
     return -1.0;
   if(i<nf2){
     for(b=0;b<nblocs;b++){val1(b)=blockF2(i,b);}
     }else{
     for(b=0;b<nblocs;b++){val1(b)=blockF3(i-nf2,b);}       
     }
     Qmat(i,i)=bjack_cov(val1,val1);
     p.increment(); 
   for(j=i+1;j<nstats;j++){
     if (Progress::check_abort() )
       return -1.0;
     if(j<nf2){
       for(b=0;b<nblocs;b++){val2(b)=blockF2(j,b);}
     }else{
       for(b=0;b<nblocs;b++){val2(b)=blockF3(j-nf2,b);}       
     }
     Qmat(i,j)=bjack_cov(val1,val2);
     Qmat(j,i)=Qmat(i,j);
     p.increment(); 
   }
 }
 for(b=0;b<nblocs;b++){val1(b)=blockF3(nf3-1,b);}
 Qmat(nstats-1,nstats-1)=bjack_cov(val1,val1);
   
 return Qmat ;
 }

//' @title compute_snpQ1onepop
//' @name compute_snpQ1onepop
//' @rdname compute_snpQ1onepop
//'
//' @description
//' Compute SNP-specific Q1 for one pop
//'
//' @param refcount Vector of nsnp counts (genotype or reads) for the reference allele
//' @param totcount Vector of nsnp total counts or read coverages
//' @param weight Numeric (w=1 for allele count data and w=poolsize/(poolsize-1) for PoolSeq data)
//'
//' @details
//' Compute SNP-specific Q1 for one pop. samples. 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q1
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ1onepop")]]
Rcpp::NumericVector compute_snpQ1onepop(Rcpp::NumericVector refcount,
                                       Rcpp::NumericVector totcount,
                                       double weight){
   int nsnps=refcount.size();
   int i ;
   double curval ;
   
   Rcpp::NumericVector Q1val(nsnps) ; // (initialise 0)
   
   Q1val.fill(NA_REAL); 
   for(i=0;i<nsnps;i++){
     if(totcount(i)>1){
       curval=(2.*refcount(i)*(totcount(i)-refcount(i)))/((totcount(i)-1.)*totcount(i)) ;
       Q1val(i)=1. - weight*curval ;
     }
   }
   return Q1val ;
 }

//' @title compute_snpQ2onepair
//' @name compute_snpQ2onepair
//' @rdname compute_snpQ2onepair
//'
//' @description
//' Compute SNP-specific Q2 for a single pair of samples
//'
//' @param refcount1 Vector of count (genotype or reads) for the reference allele in the first sample
//' @param refcount2 Vector of count (genotype or reads) for the reference allele in the second sample
//' @param totcount1 Vector of total count or read coverages in the first sample
//' @param totcount2 Vector of total count or read coverages in the second sample
//'
//' @details
//' Compute SNP-specific Q2 for a single pair of samples 
//' 
//' @return Return a vector of length nsnps with SNP-specific Q1
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".compute_snpQ2onepair")]]
Rcpp::NumericVector compute_snpQ2onepair(Rcpp::NumericVector refcount1,
                                         Rcpp::NumericVector refcount2,
                                         Rcpp::NumericVector totcount1,
                                         Rcpp::NumericVector totcount2){
   int nsnps=refcount1.size();
   int i ;
   double dumval ;
   Rcpp::NumericVector Q2val(nsnps) ; // (initialise 0)
   
   Q2val.fill(NA_REAL); 
   for(i=0;i<nsnps;i++){
     dumval=totcount1(i)*totcount2(i) ;
     if(dumval>0){
       Q2val(i)=(refcount1(i)*refcount2(i)+(totcount1(i)-refcount1(i))*(totcount2(i)-refcount2(i)))/dumval ;
     }
   }
   return Q2val ;
 }
