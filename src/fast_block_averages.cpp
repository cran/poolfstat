#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

//' @export
// [[Rcpp::export(name=".compute_Ddenom")]]
Rcpp::NumericVector compute_Ddenom( Rcpp::NumericMatrix snpQ2,
                                    Rcpp::IntegerMatrix f2idx,
                                    Rcpp::LogicalVector verbose){
  bool display_progress = is_true(any(verbose)) ;
  int nf4=f2idx.nrow();
  int nsnp=snpQ2.nrow() ;
  int i, j , idx_1 , idx_2 ;
  double dnsnp ;
  Rcpp::NumericVector Denom(nf4) ; //mean denominator
  Progress p(nf4, display_progress);
  for(i=0;i<nf4;i++){
    if (Progress::check_abort() )
      return -1.0;
    Denom(i)=0 ;
    dnsnp=0 ;
   for(j=0;j<nsnp;j++){
     idx_1=f2idx(i,0)-1 ; idx_2=f2idx(i,1)-1 ;
     if(!(NumericMatrix::is_na(snpQ2(j,idx_1)) | NumericMatrix::is_na(snpQ2(j,idx_2)))){
       dnsnp++ ;
       Denom(i)+=(1.-snpQ2(j,idx_1))*(1.-snpQ2(j,idx_2)) ;   
     }
  }
   Denom(i)/=dnsnp ; 
   p.increment(); 
  }
  return Denom ;
}

//' @export
// [[Rcpp::export(name=".compute_Q_bjmeans")]]
Rcpp::NumericVector compute_Q_bjmeans( Rcpp::NumericMatrix snpQ,
                                       Rcpp::IntegerVector snp_bj_id,
                                       Rcpp::LogicalVector verbose){
  bool display_progress = is_true(any(verbose)) ;
  int nQ=snpQ.ncol();
  int nsnp=snpQ.nrow() ;
  int nblocks=Rcpp::max(snp_bj_id) ;
  int i, j , k ;
  Rcpp::NumericMatrix bjmean(nQ,nblocks) ; //initialise a 0
  Rcpp::NumericVector bjnsnp(nblocks) ; //initialise a 0
  double sumbj,sumnsnp ;
  Progress p(nQ, display_progress);
  for(i=0;i<nQ;i++){
    if (Progress::check_abort() ){return -1.0;}
    sumbj=0 ; sumnsnp=0 ;
    //scanning snp values
    for(j=0;j<nsnp;j++){
      if(!(NumericMatrix::is_na(snpQ(j,i)))){
        k=snp_bj_id(j)-1 ;
        bjnsnp(k)++ ;
        bjmean(i,k)+=snpQ(j,i) ; 
        sumbj+=snpQ(j,i) ;
        sumnsnp++ ;
      }
    }
    //computing block values
    for(j=0;j<nblocks;j++){
      bjmean(i,j)=(sumbj-bjmean(i,j))/(sumnsnp-bjnsnp(j)) ;
      bjnsnp(j)=0; //reinitialisation
    }
    p.increment(); 
  }
  return bjmean ;
}

//' @export
// [[Rcpp::export(name=".compute_F2_bjmeans")]]
Rcpp::NumericVector compute_F2_bjmeans( Rcpp::NumericMatrix snpQ1,
                                        Rcpp::NumericMatrix snpQ2,
                                        Rcpp::IntegerMatrix q1_idx,
                                        Rcpp::IntegerVector snp_bj_id,
                                        Rcpp::LogicalVector verbose){
  bool display_progress = is_true(any(verbose)) ;
  int nF2=q1_idx.nrow();
  int nsnp=snpQ1.nrow() ; //identifcal to snpQ2.nrow() by construction
  int nblocks=Rcpp::max(snp_bj_id) ;
  int i, j , k ;
  int p1_id,p2_id ;
  Rcpp::NumericMatrix bjmean(nF2,nblocks) ; //initialise a 0
  Rcpp::NumericVector bjnsnp(nblocks) ; //initialise a 0
  double sumbj,sumnsnp,tmp_val ;
  Progress p(nF2, display_progress);
  for(i=0;i<nF2;i++){
    if (Progress::check_abort() ){return -1.0;}
    sumbj=0 ; sumnsnp=0 ;
    //scanning snp values
    for(j=0;j<nsnp;j++){
       p1_id=q1_idx(i,0)-1;
       p2_id=q1_idx(i,1)-1;      
      if(!(NumericMatrix::is_na(snpQ1(j,p1_id)) | NumericMatrix::is_na(snpQ1(j,p1_id)) | NumericMatrix::is_na(snpQ2(j,i)))){
        k=snp_bj_id(j)-1 ;
        bjnsnp(k)++ ;
        tmp_val=(snpQ1(j,p1_id)+snpQ1(j,p2_id))/2. - snpQ2(j,i) ; 
        bjmean(i,k)+=tmp_val ; 
        sumbj+=tmp_val ;
        sumnsnp++ ;
      }
    }
    //computing block values
    for(j=0;j<nblocks;j++){
      bjmean(i,j)=(sumbj-bjmean(i,j))/(sumnsnp-bjnsnp(j)) ;
      bjnsnp(j)=0; //reinitialisation
    }
    p.increment(); 
  }
  return bjmean ;
}

//' @export
// [[Rcpp::export(name=".compute_Ddenom_bjmeans")]]
Rcpp::NumericVector compute_Ddenom_bjmeans( Rcpp::NumericMatrix snpQ2,
                                            Rcpp::IntegerMatrix f2idx,
                                            Rcpp::IntegerVector snp_bj_id,
                                            Rcpp::LogicalVector verbose){
  bool display_progress = is_true(any(verbose)) ;
  int nD=f2idx.nrow();
  int nsnp=snpQ2.nrow() ; //identifcal to snpQ2.nrow() by construction
  int nblocks=Rcpp::max(snp_bj_id) ;
  int i, j , k ;
  int idx_1 , idx_2 ;
  Rcpp::NumericMatrix bjmean(nD,nblocks) ; //initialise a 0
  Rcpp::NumericVector bjnsnp(nblocks) ; //initialise a 0
  double sumbj,sumnsnp,tmp_val ;
  Progress p(nD, display_progress);
  for(i=0;i<nD;i++){
    if (Progress::check_abort() ){return -1.0;}
    sumbj=0 ; sumnsnp=0 ;
    //scanning snp values
    for(j=0;j<nsnp;j++){
      idx_1=f2idx(i,0)-1 ; idx_2=f2idx(i,1)-1 ;    
      if(!(NumericMatrix::is_na(snpQ2(j,idx_1)) | NumericMatrix::is_na(snpQ2(j,idx_2)))){
        k=snp_bj_id(j)-1 ;
        bjnsnp(k)++ ;
        tmp_val=(1.-snpQ2(j,idx_1))*(1.-snpQ2(j,idx_2)) ; 
        bjmean(i,k)+=tmp_val ; 
        sumbj+=tmp_val ;
        sumnsnp++ ;
      }
    }
    //computing block values
    for(j=0;j<nblocks;j++){
      bjmean(i,j)=(sumbj-bjmean(i,j))/(sumnsnp-bjnsnp(j)) ;
      bjnsnp(j)=0; //reinitialisation
    }
    p.increment(); 
  }
  return bjmean ;
}
