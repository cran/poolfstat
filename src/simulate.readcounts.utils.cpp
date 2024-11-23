#include <Rcpp.h>
#include <algorithm> // Required for std::shuffle
#include <random>    // Required for std::mt19937 and std::random_device
using namespace Rcpp;

//' @title simureads_poly
//' @name simureads_poly
//' @rdname simureads_poly
//'
//' @description
//' Simulate read counts from count data
//'
//' @param y_count Integer Matrix with nsnp rows and npop columns giving allele counts at the reference allele
//' @param n_count Integer Matrix with nsnp rows and npop columns giving total counts
//' @param lambda Numeric Vector of length npop giving the expected coverage of each pool
//' @param overdisp Numeric value giving overdispersion of coverages and their distribution (see details)
//' @param min_rc Integer giving the minimal read count for an allele to be considered as true allele
//' @param min_maf Float giving the MAF threshold for SNP filtering
//' @param eps Numeric value giving the sequencing error
//' @param eps_exp Numeric value giving the experimental error leading to unequal contribution of individual to the pool reads
//' @details
//'  The function implements a simulation approach similar to that described in Gautier et al. (2021). Read coverages are sampled
//'  from a distribution specified by the lambda and overdisp vectors. Note that overdisp is the same for all pop sample but 
//'  lambda (expected coverages) may vary across pool. If overdisp=1 (default in the R function), coverages are assumed Poisson distributed
//'  and the mean and variance of the coverages for the pool are both equal to the value specified in the lambda vector. If overdisp>1, coverages
//'  follows a Negative Binomial distribution with a mean equal the lamda but a variance equal to overdisp*lambda. Finally, if overdisp<1,
//'  no variation in coverage is introduced and all coverages are equal to the value specified in the lambda vector 
//'  although they may (slightly) vary in the output when eps>0 due to the removal of error reads.
//'  The eps parameter control sequencing error rate. Sequencing errors are modeled following Gautier et al. (2021) i.e. read counts for the four
//'  possible bases are sampled from a multinomial distribution Multinom(c,\{f*(1-eps)+(1-f)*eps/3;f*eps/3+(1-f)*(1-eps),eps/3,eps/3\}) 
//'  where c is the read coverage and f the reference allele frequencies (obtained from the count data).
//'  Experimental error eps_exp control the contribution of individual (assumed diploid) to the pools following the model described 
//'  in Gautier et al. (2013).  The parameter eps_exp corresponds to the coefficient of variation of the individual contributions
//'  When eps_exp tends toward 0, all individuals contribute equally to the pool and there is no experimental error. For example, 
//'  with 10 individuals, eps_exp=0.5 correspond to a situation where 5 individuals contribute 2.8x more reads than the five others.
//'  Note that the number of (diploid) individuals for each SNP and pop. sample is deduced from the input total count 
//'  (it may thus differ over SNP when the total counts are not the same). 
//'  
//' @return Return an Integer matrix with nsnp rows and 2*npop columns (1:npop=ref allele readcount; (npop+1):2*npop=coverage)  
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".simureads_poly")]]
IntegerMatrix simureads_poly(IntegerMatrix y_count,
                             IntegerMatrix n_count,
                             NumericVector lambda,
                             double overdisp,
                             int min_rc,
                             double min_maf,
                             double eps,
                             double eps_exp){
  RNGScope scope; // Initialize RNG within this scope
  int npos = y_count.nrow() ;
  int npop = y_count.ncol() ;
  int nall=4,ndip,cur_cov ;
  double ref_freq,p_nbinom,rho ;
  Rcpp::IntegerVector allele_idx=seq(0,nall-1) ;
  NumericVector probs_multi,ind_contrib,r_nbinom(npop);
  IntegerMatrix readcounts(npos,2*npop);  //nread1 pour pop 1 a npop; coverage pour pop 1 a npop
  IntegerMatrix cur_counts(npop,nall) ;
  bool nbinom=FALSE,nosampling=FALSE,exp_error=FALSE,seq_error=FALSE;
  if(eps>0){seq_error=TRUE;}
  if(overdisp<1){
    nosampling=TRUE;
  }else{
   if(overdisp>1){
     nbinom=TRUE; 
     p_nbinom=1./overdisp ; r_nbinom=lambda/(overdisp-1.); 
   }
   }
  if(eps_exp>1e-3){exp_error=TRUE;}
  
  for(int i=0;i<npos;i++){
  //  Rcpp::Rcout << i+1 << "\n" ;
    IntegerVector tmp_all_count(nall) ;
    //position monomorphe
    for(int k=0;k<npop;k++){
      if(n_count(i,k)==0){
        continue; //not necessary to evaluate further!
      }
    // sampling coverage
      if(nosampling){
        cur_cov=lambda(k) ;
      }else{
        if(nbinom){
          cur_cov = rnbinom(1,r_nbinom(k),p_nbinom)[0] ; 
        }else{
          cur_cov = rpois(1,lambda(k))[0] ; 
      }
      }
      if(cur_cov==0){
        continue; //not necessary to evaluate further!
      }
      // 
      IntegerVector out(nall) ;
      if(exp_error){
        ndip=floor(double(n_count(i,k))/2.) ;
        int tot_refcount=0; //count of ref allele reads over all individuals
        if(ndip>0){
          if(ndip==1){
            if(y_count(i,k)==1){
              tot_refcount=rbinom(1,cur_cov,0.5)[0];
            }else{
              if(y_count(i,k)==2){tot_refcount=cur_cov;} 
            }            
          }else{
           IntegerVector out_indcov(ndip) ;
           rho=(double(ndip) - 1.)/(eps_exp*eps_exp) - 1. ; 
           //sample from a dirichlet
           NumericVector ind_contrib=Rcpp::rgamma(ndip,rho/double(ndip),1.) ;
           ind_contrib = ind_contrib/sum(ind_contrib);
           rmultinom(cur_cov, ind_contrib.begin(), ndip, out_indcov.begin());
           //vector of individual counts by shuffling pool count
           IntegerVector indcount(n_count(i,k));
           std::fill(indcount.begin(), indcount.begin() + y_count(i,k), 1);
           std::fill(indcount.begin() + y_count(i,k), indcount.end(), 0);
           // Shuffle the vector using std::shuffle (better than random_shuffle deprecated)
           std::shuffle(indcount.begin(), indcount.end(), std::mt19937(std::random_device()()));
           //count the number of ref allele reads over all diploid individuals
           for(int id=0;id<ndip;id++){
             if(out_indcov(id)==0){
               continue;
             }
             int id_allidx=2*id;
             if(indcount(id_allidx)!=indcount(id_allidx+1)){
               tot_refcount+=rbinom(1,out_indcov(id),0.5)[0] ;
              }else{
               if(indcount(id_allidx)==1){tot_refcount+=out_indcov(id) ;}
           }
          }
          }
          //define the vector of pop read for all alleles
          if(seq_error){
           ref_freq=double(tot_refcount)/double(cur_cov) ;
           probs_multi={(1.-eps)*ref_freq+(1.-ref_freq)*eps/3.,(1.-eps)*(1.-ref_freq)+ref_freq*eps/3.,eps/3.,eps/3.} ;
           rmultinom(cur_cov, probs_multi.begin(), nall, out.begin());            
           }else{out(0)=tot_refcount ; out(1)=cur_cov-tot_refcount;}
        }else{//ndip==0
          continue;}
      }else{//no experimental error
        ref_freq=double(y_count(i,k))/double(n_count(i,k)) ;
        if(seq_error){
         probs_multi={(1.-eps)*ref_freq+(1.-ref_freq)*eps/3.,(1.-eps)*(1.-ref_freq)+ref_freq*eps/3.,eps/3.,eps/3.} ;
         rmultinom(cur_cov, probs_multi.begin(), nall, out.begin());
         }else{
           out(0)=rbinom(1,double(cur_cov),ref_freq)[0];//R::rbinom(double(cur_cov),ref_freq);
           out(1)=cur_cov-out(0);
         }
      }
      cur_counts(k,_)=out ;
      for(int l=0;l<nall;l++){
          tmp_all_count[l]+=cur_counts(k,l) ;
        }
    }
    
    if(seq_error){
     std::sort(allele_idx.begin(), allele_idx.end(), [&](int a, int b){
     return tmp_all_count[a] > tmp_all_count[b];});//on ordonne les index
     int refall1=allele_idx(0) ;
     int refall2=allele_idx(1) ;
  //      Rcpp::Rcout << i+1 << " " << tmp_all_count <<"\n" ;
       if(tmp_all_count[allele_idx(2)]<=min_rc && //le 4 l'est forcement!
          tmp_all_count[refall2]>=min_rc){
       double maf= double (tmp_all_count[refall2])/(double (tmp_all_count[refall2])+ double (tmp_all_count[refall1])) ; 
   //   Rcpp::Rcout << i+1 << " " << tmp_all_count[refall1] << " " << tmp_all_count[refall2] << " " << maf<<"\n" ;  
       if(maf>min_maf){
         for(int k=0;k<npop;k++){readcounts(i,npop+k)=cur_counts(k,refall1) + cur_counts(k,refall2) ;}
         if(refall1==0 || refall2==0){//the original ref allele is the major one or the minor one (we keep the same coding)
          for(int k=0;k<npop;k++){readcounts(i,k)=cur_counts(k,0) ;}
         }else{//the original ref allele is no more represented (may be rare!): in this case the ref allele is the minor one
           for(int k=0;k<npop;k++){readcounts(i,k)=cur_counts(k,refall2) ;}           
         }  
      }
      }
    }else{
      double maf=double(tmp_all_count(0))/(double(tmp_all_count(0))+double(tmp_all_count(1)));
      if(maf>0.5){maf=1-0.5;}
      if(tmp_all_count(0)>=min_rc && tmp_all_count(1)>=min_rc && maf>min_maf){
        for(int k=0;k<npop;k++){
          readcounts(i,k)=cur_counts(k,0) ;
          readcounts(i,npop+k)=cur_counts(k,0) + cur_counts(k,1) ;          
        }
      }
    }
    }
  return readcounts;
}

//' @title simureads_mono
//' @name simureads_mono
//' @rdname simureads_mono
//'
//' @description
//' Simulate read counts for monomorphic position when there is sequencing error
//'
//' @param npos Integer giving the number of positions (close to genome size)
//' @param npop Integer giving the number of population samples
//' @param lambda Numeric Vector of length npop giving the expected coverage of each pool
//' @param overdisp Numeric value giving overdispersion of coverages and their distribution (see details)
//' @param min_rc Integer giving the minimal read count for an allele to be considered as true allele
//' @param min_maf Float giving the MAF threshold for SNP filtering
//' @param eps Numeric value giving the sequencing error
//' @details
//' The function implements a simulation approach similar to that described in Gautier et al. (2021). Read coverages are sampled
//' from a distribution specified by the lambda and overdisp vectors. Note that overdisp is the same for all pop sample but 
//' lambda (expected coverages) may vary across pool. If overdisp=1 (default in the R function), coverages are assumed Poisson distributed
//' and the mean and variance of the coverages for the pool are both equal to the value specified in the lambda vector. If overdisp>1, coverages
//' follows a Negative Binomial distribution with a mean equal the lamda but a variance equal to overdisp*lambda. Finally, if overdisp<1,
//' no variation in coverage is introduced and all coverages are equal to the value specified in the lambda vector 
//' although they may (slightly) vary in the output when eps>0 due to the removal of error reads.
//' The eps parameter control sequencing error rate. Sequencing errors are modeled following Gautier et al. (2021) i.e. read counts for the four
//' possible bases are sampled from a multinomial distribution Multinom(c,\{1-eps;eps/3,eps/3,eps/3\}) 
//' where c is the read coverage. Only bi-allelic SNPs (after considering min_rc) satisfying with MAF>min_maf are included in the output.
//'  
//' @return Return an Integer matrix with nsnp rows and 2*npop columns (1:npop=ref allele readcount; (npop+1):2*npop=coverage)  
//' 
//' @examples
//' #
//' @export
// [[Rcpp::export(name=".simureads_mono")]]
IntegerMatrix simureads_mono(int npos,
                             int npop,
                             NumericVector lambda,
                             double overdisp,
                             int min_rc, //nbre minimal pour mas etre exclu e.g., 2
                             double min_maf,
                             double eps){
  RNGScope scope; // Initialize RNG within this scope
  int cur_cov,nsnp=0,dum ;
  double ref_freq,p_nbinom,r_nbinom,lambda_cov ;
  Rcpp::IntegerVector nerr,outcount,outcoverage,popcov(npop) ;
  bool nbinom=FALSE,nosampling=FALSE;
  if(overdisp<1){
    nosampling=TRUE;
  }else{
   if(overdisp>1){
     nbinom=TRUE; 
     p_nbinom=1./overdisp ; 
     //Sum of n neg binomial of same p NB(r,p) follows NB(nr,p)
     r_nbinom=sum(lambda)/(overdisp-1.); 
   }else{
     //Sum of n Poisson(l_i) follows Poisson(Sum(l_i))
     lambda_cov=sum(lambda);
   }
   }
  //distribution des couvertures par pops (evite de perdre du temps a echantilloner des couvertures + permet de mieux gere la memoire)
  int ncov_values;
  IntegerVector cov_values,cov_nsamples;
  if(nosampling){
    ncov_values=1 ;
    cov_values.push_back(lambda_cov);cov_nsamples.push_back(npos);
  }else{
    int mincov,maxcov;
    if(nbinom){
      mincov=R::qnbinom(1e-5,r_nbinom,p_nbinom,true,false) ;
      maxcov=R::qnbinom(1e-5,r_nbinom,p_nbinom,false,false) ;      
    }else{
      mincov=R::qpois(1e-5,lambda_cov,true,false) ;      
      maxcov=R::qpois(1e-5,lambda_cov,false,false) ;
    }
 //   Rcpp::Rcout << "mincov: " << mincov << " maxcov: " << maxcov << " lda: " << lambda_cov <<" POS\n" ;
    cov_values=seq(mincov,maxcov) ;
    ncov_values=cov_values.length() ;
    for(int k=0;k<ncov_values;k++){
      if(nbinom){
        cov_nsamples.push_back(round(R::dnbinom(cov_values[k],r_nbinom,p_nbinom,false)*double(npos)));        
    }else{
        cov_nsamples.push_back(round(R::dpois(cov_values[k],lambda_cov,false)*double(npos)));
    }
  }
     npos=sum(cov_nsamples);
    }
  
 //  Rcpp::Rcout << npos << " " << cov_nsamples.length() << " " <<cov_values.length() <<" POS\n" ;
  
  //echantillonnage des errors
  int countalt,naltok,tmp_readtot;
  IntegerVector nread_tot,nread_error,nalt(3);
  IntegerVector nsamp_reads_err,samp_read_err,samp_read_tot; //used for count algo  
  NumericVector prob_err(3);
  prob_err.fill(1./3.) ;
  for(int c=0;c<ncov_values;c++){
    if(cov_values(c)==0 || cov_nsamples(c)==0){continue;}
 //   Rcpp::Rcout << c << " out of " << ncov_values << " : "  <<cov_values(c)<< " " << cov_nsamples(c) <<"\n" ;
    if(cov_nsamples(c)<100){//direct sampling (may be inefficient for high value if lot of case OK)
     nerr=rbinom(cov_nsamples(c),cov_values(c),eps);
     dum=floor(double(cov_values(c))*min_maf) ;
     for(int i=0;i<cov_nsamples(c);i++){
      if(nerr(i)==0){continue;}
      if(nerr(i)<min_rc){continue;}     
      if(nerr(i)>dum){
       rmultinom(nerr(i), prob_err.begin(), 3, nalt.begin());
       countalt=max(nalt);
       if(countalt>=min_rc &&  (double(countalt)/double(tmp_readtot))>min_maf){//only bi-allelic according to conditions on min_rc to retain alleles as ture alleles
       naltok=0;
 //     Rcpp::Rcout << c << " "  <<cov_values(c)<< " " << cov_nsamples(c) <<" "<<countalt<<"\n" ;
       for(int k=0;k<3;k++){
         if(nalt(k)>=min_rc){naltok++;}
       }
       tmp_readtot=cov_values(c)-sum(nalt)+countalt;
       if(naltok==1){ 
         nread_tot.push_back(tmp_readtot);
         nread_error.push_back(countalt);
         nsnp++;
       }
     }
    }
     }
  }else{//counts from expectations
    int nerr_max=R::qbinom(1./cov_nsamples(c),cov_values(c),eps,false,false) ;
    int min_e=floor(double(cov_values(c))*min_maf) ;
    if(min_e<min_rc){min_e=min_rc ;}
    if(nerr_max>=min_rc){
     for(int e=nerr_max;e>=min_rc;e--){
       int nsamp=ceil(R::dbinom(e,cov_values(c),eps,false) *cov_nsamples(c)) ;
       if(nsamp>0){
       for(int ei=e;ei>=min_e;ei--){//nber of alt all (for conf OK)
        int nothers=e-ei ; //the two latest allele
         if(nothers>=2*min_rc){continue;} //one or the 2 other alleles have count>min_rc => tri-allelic
         for(int ej=0;ej<=nothers;ej++){
           if(ej<min_rc){
           for(int ek=0;ek<=nothers;ek++){
             if(ek<min_rc && (ej+ek)==nothers){// conf OK
              double multinom_pmf = lgamma(double(e) + 1) + log(1./3.)*double(e) - 
                (lgamma(double(ei)+1.) + lgamma(double(ej)+1.) + lgamma(double(ek)+1.));
              multinom_pmf=exp(multinom_pmf);
              int noccurence=ceil(double(nsamp)*multinom_pmf*3.); //*3 because 3 possible way of chosing alt allele
 //    Rcpp::Rcout << nerr_max << " " << ei << " , "  << ej << " , "  << ek << " , " << noccurence << " " << nsamp << " " << multinom_pmf << "\n" ;
               if(noccurence>0){//new vectors not using push_back => too slow if large vectors
                 nsamp_reads_err.push_back(noccurence);
                 samp_read_err.push_back(ei) ;
                 samp_read_tot.push_back(cov_values(c)-ej-ek); 
               }
              }
            }
          }
           }
        }}
      }
    }
  }
  }
  //concatenation counts from sampling and from expectations
  int size2=sum(nsamp_reads_err) ;
  if(size2>0){
    nsnp+=size2 ;
    int ncateg=nsamp_reads_err.size() ;
    int size1=nread_tot.size() ;
    IntegerVector new_vec1(size1+size2),new_vec2(size1+size2); //new_read_tot et nread_error
    for(int k=0;k<size1;k++){//from sampling
      new_vec1(k)=nread_tot(k) ; new_vec2(k)=nread_error(k) ;
    }
    int dum=0 ;
    for(int k=0;k<ncateg;k++){
      for(int l=0;l<nsamp_reads_err(k);l++){
        new_vec1(dum+size1)=samp_read_tot(k) ;
        new_vec2(dum+size1)=samp_read_err(k) ; 
        dum++ ;
      }
    }
    nread_tot=new_vec1 ;
    nread_error=new_vec2 ;
  }
  
//  Rcpp::Rcout << nsnp << " new SNPs\n" ;
  
  //create counts
  IntegerMatrix readcounts(nsnp,2*npop);
  if(nsnp>0){
    NumericVector probs_multi(npop,1./(double (npop))) ; //pour position des pop (refractionner couverture)
    IntegerVector cov_pop(npop),err_indices,pos_err;
    int err_idx,id_start,id_end;
    for(int i=0;i<nsnp;i++){
      if(nread_tot(i)<1){continue;}
      rmultinom(nread_tot(i), probs_multi.begin(), npop, cov_pop.begin());
      err_indices = seq(0,nread_tot(i)-1);  // Create a sequence from 0 to n-1
      pos_err = Rcpp::sample(err_indices, nread_error(i), false);
      std::sort(pos_err.begin(), pos_err.end());
  //    Rcpp::Rcout << pos_err<< "\n" ; Rcpp::Rcout << cov_pop<< "\n" ;      
      err_idx=0;id_start=0;
      for(int k=0;k<npop;k++){
        if(cov_pop(k)<1){continue;}
        id_end=id_start+cov_pop(k)-1;
        readcounts(i,npop+k)=cov_pop(k) ;
        readcounts(i,k)=cov_pop(k) ;
        while(err_idx<nread_error(i) && pos_err(err_idx)<=id_end){
          readcounts(i,k)--;
          err_idx++;
        }
        id_start=id_end+1;
      }
      }
    }
  return readcounts;  
}

