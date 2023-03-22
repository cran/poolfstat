#include <Rcpp.h>
#include <string>

using namespace Rcpp;
//using namespace string;

//' @title scan_allele_info
//' @name scan_allele_info
//' @rdname scan_allele_info
//'
//' @description
//' Scan allele information in ALT field of a vcf
//'
//' @param allele_info a character string vector (ALT field of the vcf)
//'
//' @details
//' Scan allele information in ALT field of a vcf to identify the number of alleles and if there is indels
//'
//' @return Return a vector with two elements consisting i) the number of alleles (1+number of comma)
//' and ii) 0 or 1 if an indel is detected 
//' 
//' @examples
//' .scan_allele_info(c("A,C","T","AAT"))
//' 
//' @export
// [[Rcpp::export(name=".scan_allele_info")]]
Rcpp::IntegerMatrix scan_allele_info( Rcpp::StringVector allele_info) {
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


//' @title extract_vscan_counts
//' @name extract_vscan_counts
//' @rdname extract_vscan_counts
//'
//' @description
//' Extract VarScan counts
//'
//' @param vcf_data a matrix of String containing count information in VarScan format
//' @param ad_idx the index of the FORMAT AD field  
//' @param rd_idx the index of the FORMAT RD field 
//'
//' @details Extract VarScan counts and return read counts for the reference and alternate allele.
//' For VarScan generated vcf, SNPs with more than one alternate allele are discarded 
//' (because only a single count is then reported in the AD fields) making the min.rc unavailable (of vcf2pooldata).
//' The VarScan --min-reads2 option might replace to some extent the min.rc functionality although 
//' SNP where the two major alleles in the Pool-Seq data are different from the reference allele 
//' (e.g., expected to be more frequent when using a distantly related reference genome for mapping) 
//' will be disregarded.
//' @return A numeric matrix of read count with nsnp rows and 2*npools columns.
//' The first npools columns consist of read count for the reference allele (RD),
//' columns npools+1 to 2*npools consist of read coverage (RD+AD)
//' @examples
//' .extract_vscan_counts(rbind(c("0/0:0:20","1/1:18:1"),c("0/1:12:15","1/1:27:2")),3,2)
//' 
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

//' @title extract_nonvscan_counts
//' @name extract_nonvscan_counts
//' @rdname extract_nonvscan_counts
//'
//' @description
//' Extract counts from vcf produced by other caller than VarScan (e.g., bcftools, FreeBayes, GATK)
//'
//' @param vcf_data a matrix of String containing count information
//' @param nb_all a vector containing the number of alleles for the different markers
//' @param ad_idx the index of the FORMAT AD field  
//' @param min_rc Minimal allowed read count per base (same as min.rc option in \code{\link{vcf2pooldata}}) 
//'
//' @details Extract VarScan counts and return read counts for the reference and alternate allele
//' @return A numeric matrix of read count with nsnp rows and 2*npools+6 columns.
//' The first npools columns consist of read count for the reference allele,
//' columns npools+1 to 2*npools consist of read coverage. The last 6 columns correspond to 
//' the index of the two most frequent alleles (idx_all1 and idx_all2) and their count (cnt_all1 and cnt_all2);
//' the min_rc filtering criterion and count of variant (cnt_bases) other than two first most frequent. The min_rc crit is
//' set to -1 for polymorphisms with more than 2 alleles and with the third most frequent alleles having 
//' more than min_rc count 
//' @examples
//' .extract_nonvscan_counts(rbind(c("0/0:20,0","1/1:1,18"),c("0/2:12,1,15","1/1:27,1,0")),c(2,3),2,0)
//' .extract_nonvscan_counts(rbind(c("0/0:20,0","1/1:1,18"),c("0/2:12,1,15","1/1:27,1,0")),c(2,3),2,2)
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
                                                      // min_rc crit is set to -1 for polymorphisms with more than 2 alleles with the third most frequent alleles having more than min_rc count 
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



//' @title extract_allele_names
//' @name extract_allele_names
//' @rdname extract_allele_names
//'
//' @description
//' Extract the alleles from the REF and ALT fields
//'
//' @param allele_info a character string vector (concatenated REF and ALT field of the vcf)
//' @param allele_idx Matrix with indexes of the two alleles of interest for the different markers
//'
//' @details
//' Extract the alleles from the REF and ALT fields
//' 
//' @return Return a matrix with the two alleles after parsing the alleles info
//' 
//' @examples
//' .extract_allele_names(c("A,C","A,C,T"),rbind(c(1,2),c(1,3)))
//' 
//' @export
// [[Rcpp::export(name=".extract_allele_names")]]
Rcpp::StringMatrix extract_allele_names( Rcpp::StringVector allele_info,
                                         Rcpp::IntegerMatrix allele_idx ) {
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


//' @title find_indelneighbor_idx
//' @name find_indelneighbor_idx
//' @rdname find_indelneighbor_idx
//'
//' @description
//' Search for the closest indels of the markers
//'
//' @param contig a character string vector corresponding to the CHR field value of the vcf for the markers
//' @param position an integer vector corresponding to the POSITION value for the markers 
//' @param indels_idx vector of (0-indexed) indices of indels
//' @param min_dist same as min.dist.from.indels option in \code{\link{vcf2pooldata}}
//' @param indels_size size of the indels (associated to indels_idx)
//'
//' @details
//' Identify if the SNPs are close to an indel
//' 
//' @return Return a vector consisting of 1 (if the marker is close to an indel) or 0 (if not)
//' 
//' @examples
//' .find_indelneighbor_idx(c("chr1","chr1","chr1"),c(1000,1004,1020),1,5,2)
//' 
//' @export
// [[Rcpp::export(name=".find_indelneighbor_idx")]]
Rcpp::IntegerVector find_indelneighbor_idx( Rcpp::StringVector contig,
                                            Rcpp::IntegerVector position,
                                            Rcpp::IntegerVector indels_idx,
                                            int min_dist,
                                            Rcpp::IntegerVector indels_size) {//return 1 if close to indel, 0 else
  int idx_j=0 ;
  int npos=position.size();
  int nindels=indels_idx.size();
  std::string cur_ctg ;
  int cur_pos_thr;
  Rcpp::IntegerVector return_vector(npos);
  for(int i=0;i<nindels;i++){
    return_vector(indels_idx(i))=1 ;
    cur_ctg=contig(indels_idx(i));
    //to the left
    idx_j=indels_idx(i)-1;
    cur_pos_thr=position(indels_idx(i))-min_dist;
    while(idx_j>=0 && contig[idx_j]==cur_ctg && position[idx_j]>=cur_pos_thr){
      return_vector(idx_j)=1;
      idx_j--;
    }
    //to the right
    idx_j=indels_idx(i)+1;
    cur_pos_thr=position(indels_idx(i))+min_dist+indels_size(i)-1 ;
    while(idx_j<=npos && contig[idx_j]==cur_ctg && position[idx_j]<=cur_pos_thr){
      return_vector(idx_j)=1;
      idx_j++;
    }
  }
  return return_vector;
}
