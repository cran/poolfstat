#include <Rcpp.h>
#include <string>

using namespace Rcpp;
//using namespace string;

//' @export
// [[Rcpp::export(name=".scan_allele_info")]]
Rcpp::IntegerMatrix scan_allele_info( Rcpp::StringVector allele_info) {//return a two vector with number of alleles (1+number of,) and 0 or 1 if indel
  int pos_j=0 ;
  int npos=allele_info.size();
  std::string istring ;
  Rcpp::IntegerMatrix return_matrix(npos, 2);// count ref_allele 1:npools; coverage npools+1:2npools
  for(int i=0;i<npos;i++){
    istring = allele_info(i) ;
    return_matrix(i,0)=1 ; //Initialisation since there is at least the reference allele (note the ref allele may be not covered (this is treated in downstream analyses for non varscan file) 
    pos_j=-1 ; // to allow proper test for indels for the first element: e.g., if ALT=AG, first comma at pos j=2..
    for(int j=0;j<istring.size();j++){
      if(istring[j] == ','){
        return_matrix(i,0)++;
        if(j-pos_j>2){return_matrix(i,1)=1;}
        pos_j=j ;
      }
      //check last allele if indel
      if(j-pos_j>1){return_matrix(i,1)=1;}
    }
  }
  return return_matrix;
}

//' @export
// [[Rcpp::export(name=".extract_vscan_counts")]]
Rcpp::NumericMatrix extract_vscan_counts( Rcpp::StringMatrix vcf_data,
                                          int ad_idx,
                                          int rd_idx) {
  int i = 0;
  int j = 0;
  int pos_j=1 ;
  std::string istring,rd_string,ad_string ;
  int max_idx=std::max(ad_idx,rd_idx) ;
  int npos=vcf_data.nrow();
  int npools=vcf_data.ncol()  ;
//   Rcpp::Rcout << vcf_data.nrow() << " " << vcf_data.ncol() << "\n";
   Rcpp::NumericMatrix return_matrix(npos, 2*npools);// count ref_allele 1:npools; coverage npools+1:2npools
   for(int m=0;m<npos;m++){
      for(int n=0;n<npools;n++){      
       istring = vcf_data(m,n) ;
       istring.push_back(':') ; // to ensure parsing until the end (i.e., if only the field finishes by RD or AD) 
       pos_j=-1 ; j=0 ; i=1 ; //(field since should never start by ':')
       while(i < istring.size() && j!=max_idx){
         if(istring[i] == ':'){
         if(j==(ad_idx-1)){
          ad_string=istring.substr( pos_j + 1, i - pos_j - 1 ) ;
          }
         if(j==(rd_idx-1)){
          rd_string=istring.substr( pos_j + 1, i - pos_j - 1 ) ;
        }
      j++ ;  pos_j=i ;      
    }
    i++ ;
   }
  if(ad_string[0] != '.' && rd_string[0] != '.' && j==max_idx){//the last to ensure proper field
     return_matrix(m,n)=atof(rd_string.c_str()) ;
     return_matrix(m,npools+n)=atof(ad_string.c_str()) + return_matrix(m,n) ;     
  }
 //   Rcpp::Rcout << istring << " " << ad_string << " " << rd_string  << "\n" ; //<< pos1 << " " << pos2 ;
        }
   }
  return return_matrix;
}


//' @export
// [[Rcpp::export(name=".extract_nonvscan_counts")]]
Rcpp::NumericMatrix extract_nonvscan_counts( Rcpp::StringMatrix vcf_data,
                                             Rcpp::IntegerVector nb_all,
                                             int ad_idx,
                                             int min_rc) {
  int i, j, pos_j ;
  int cntall1, cntall2 ;
  std::string istring,ad_string ;
  int npos=vcf_data.nrow();
  int npools=vcf_data.ncol()  ;
  int nb_all_max=Rcpp::max(nb_all) ;
  //   Rcpp::Rcout << vcf_data.nrow() << " " << vcf_data.ncol() << "\n";
  Rcpp::NumericMatrix return_matrix(npos, 2*npools+6);// count ref_allele 1:npools; coverage npools+1:2npools; 
                                                      //idx_all1, idx_all2, cnt_all1, cnt_all2, min_rc crit, cnt_bases other than 1 and 2
                                                      // min_rc crit is set to -1 for polymorphisms with more than 2 alleles with the third most frequent alleles haveing more than min_rc count 
  Rcpp::NumericMatrix tmp_cnt_alleles(npools+1,nb_all_max) ; //Only useful if nb_all>2 (last rows for the sums over pools)
  Rcpp::IntegerVector allele_idx=seq(0,nb_all_max-1) ;
  for(int m=0;m<npos;m++){
 //   Rcpp::Rcout << m << " " << nb_all(m) << "\n";
    if(nb_all(m)>2){tmp_cnt_alleles=0*tmp_cnt_alleles ;}
    cntall1=0 ; cntall2=0 ;
    ad_string=istring[0] ;
    for(int n=0;n<npools;n++){      
      istring = vcf_data(m,n) ;
      istring.push_back(':') ; //in case ad field is at the end
      pos_j=-1 ; j=0 ; i=1 ; //(field since should never start by ':')
      while(i < istring.size() && j!=ad_idx){
  //      Rcpp::Rcout <<istring << " " << istring.size() << " " << i << " " << j <<  " " << pos_j << " "<< ad_idx << " " << ad_string <<"\n";
        if(istring[i] == ':'){
          if(j==(ad_idx-1)){
   //         Rcpp::Rcout << pos_j + 1 << " " << i-pos_j-1 << " " << istring[pos_j+1] << " " << istring.substr( pos_j + 1, i - pos_j - 1 ) << "\n";
            ad_string = istring.substr( pos_j + 1, i - pos_j - 1 ) ;
          }
          j++ ; pos_j=i ;      
        }
        i++ ;
      }
  //    if(nb_all(m)==3){ Rcpp::Rcout << m << " " << n << " " <<ad_string << " " << nb_all(m) <<  "\n" ;}
      if(ad_string[0] != '.' && j==ad_idx){//parse counts
        i=1 ; //(field since should never start by ',')
        if(nb_all(m)==2){
          while(i < ad_string.size() && ad_string[i] != ','){i++;}
          istring=ad_string.substr( 0, i ) ;
          return_matrix(m,n)=atof(istring.c_str()) ;//cnt for all1
          istring=ad_string.substr( i+1, ad_string.size() - i - 1 ) ;
          return_matrix(m,npools+n)=atof(istring.c_str()) ;//cnt for all2
          cntall1+=return_matrix(m,n) ; cntall2+=return_matrix(m,npools+n) ;
          return_matrix(m,npools+n)+=return_matrix(m,n) ; //transform into coverage
        }else{//we need to store the counts to identify the two most frequent alleles (over all pools)
          pos_j=-1 ; j=0 ; 
          ad_string.push_back(',') ; //trick to ensure parsing until the end
           while(i < ad_string.size()){
            if(ad_string[i] == ','){
              istring=ad_string.substr( pos_j + 1, i - pos_j - 1 ) ;
              tmp_cnt_alleles(n,j)=atof(istring.c_str()) ;
              tmp_cnt_alleles(npools,j)+=tmp_cnt_alleles(n,j) ;
     //         Rcpp::Rcout << m << " " << n << " " << ad_string << " " << istring << " " << tmp_cnt_alleles(n,j)  <<  "\n" ;
              j++ ;  pos_j=i ;
            }
            i++;
            } 
        }
        }
    }
    if(nb_all(m)==2){
        return_matrix(m,2*npools)=1 ; return_matrix(m,2*npools+1)=2 ;      
        return_matrix(m,2*npools+2)=cntall1 ; return_matrix(m,2*npools+3)=cntall2 ;        
    }else{
      std::sort(allele_idx.begin(), allele_idx.end(), [&](int a, int b){return tmp_cnt_alleles(npools,a) > tmp_cnt_alleles(npools,b);});//on ordonne les index
      i=allele_idx(0) ; j=allele_idx(1) ;
      return_matrix(m,2*npools)=i+1 ; return_matrix(m,2*npools+1)=j+1 ;
      return_matrix(m,2*npools+2)=tmp_cnt_alleles(npools,i) ; 
      return_matrix(m,2*npools+3)=tmp_cnt_alleles(npools,j) ;
      if(tmp_cnt_alleles(npools,allele_idx(2))>min_rc){
        return_matrix(m,2*npools+4)=-1 ;
      }
      return_matrix(m,2*npools+5)=sum(tmp_cnt_alleles(npools,_)) - return_matrix(m,2*npools+2) - return_matrix(m,2*npools+3) ;  
      // Rcpp::Rcout <<tmp_cnt_alleles(npools,i)  << " " << tmp_cnt_alleles(npools,j) << " " << return_matrix(m,2*npools+5) << " " << i << " " << j <<  "\n" ;
             for(int n=0;n<npools;n++){ //fill the matrix with the two major alleles 
         return_matrix(m,n)=tmp_cnt_alleles(n,i) ;
         return_matrix(m,npools+n)=tmp_cnt_alleles(n,i) + tmp_cnt_alleles(n,j);
      }
    }
      //   Rcpp::Rcout << istring << " " << ad_string << " " << rd_string  << "\n" ; //<< pos1 << " " << pos2 ;
    }
  
  return return_matrix;
}

//' @export
// [[Rcpp::export(name=".extract_allele_names")]]
Rcpp::StringMatrix extract_allele_names( Rcpp::StringVector allele_info,
                                         Rcpp::IntegerMatrix allele_idx ) {//return a two-dim matrix with the two alleles after parsing the alleles info
  int npos=allele_info.size();
  std::string istring ;
  Rcpp::StringMatrix return_matrix(npos, 2);// allele name for first and second idx
  for(int i=0;i<npos;i++){
    istring = allele_info(i) ;
    int pos_j=-1 ;
    int j=0 ; 
    int k=0 ;
    Rcpp::StringVector dum_str ;
    istring.push_back(',') ; //trick to ensure parsing until the end
    while(k < istring.size()){
      if(istring[k] == ','){
        dum_str.push_back(istring.substr( pos_j + 1, k - pos_j - 1 )) ;
        //         Rcpp::Rcout << m << " " << n << " " << ad_string << " " << istring << " " << tmp_cnt_alleles(n,j)  <<  "\n" ;
        j++ ;  pos_j=k ;
      }
      k++;
    }
    return_matrix(i,0)=dum_str(allele_idx(i,0)-1) ;
    return_matrix(i,1)=dum_str(allele_idx(i,1)-1) ;
    }
  return return_matrix;
}
