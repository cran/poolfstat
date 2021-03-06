// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// scan_allele_info
Rcpp::IntegerMatrix scan_allele_info(Rcpp::StringVector allele_info);
RcppExport SEXP _poolfstat_scan_allele_info(SEXP allele_infoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type allele_info(allele_infoSEXP);
    rcpp_result_gen = Rcpp::wrap(scan_allele_info(allele_info));
    return rcpp_result_gen;
END_RCPP
}
// extract_vscan_counts
Rcpp::NumericMatrix extract_vscan_counts(Rcpp::StringMatrix vcf_data, int ad_idx, int rd_idx);
RcppExport SEXP _poolfstat_extract_vscan_counts(SEXP vcf_dataSEXP, SEXP ad_idxSEXP, SEXP rd_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix >::type vcf_data(vcf_dataSEXP);
    Rcpp::traits::input_parameter< int >::type ad_idx(ad_idxSEXP);
    Rcpp::traits::input_parameter< int >::type rd_idx(rd_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_vscan_counts(vcf_data, ad_idx, rd_idx));
    return rcpp_result_gen;
END_RCPP
}
// extract_nonvscan_counts
Rcpp::NumericMatrix extract_nonvscan_counts(Rcpp::StringMatrix vcf_data, Rcpp::IntegerVector nb_all, int ad_idx, int min_rc);
RcppExport SEXP _poolfstat_extract_nonvscan_counts(SEXP vcf_dataSEXP, SEXP nb_allSEXP, SEXP ad_idxSEXP, SEXP min_rcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix >::type vcf_data(vcf_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nb_all(nb_allSEXP);
    Rcpp::traits::input_parameter< int >::type ad_idx(ad_idxSEXP);
    Rcpp::traits::input_parameter< int >::type min_rc(min_rcSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_nonvscan_counts(vcf_data, nb_all, ad_idx, min_rc));
    return rcpp_result_gen;
END_RCPP
}
// extract_allele_names
Rcpp::StringMatrix extract_allele_names(Rcpp::StringVector allele_info, Rcpp::IntegerMatrix allele_idx);
RcppExport SEXP _poolfstat_extract_allele_names(SEXP allele_infoSEXP, SEXP allele_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type allele_info(allele_infoSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type allele_idx(allele_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_allele_names(allele_info, allele_idx));
    return rcpp_result_gen;
END_RCPP
}
// compute_Ddenom
Rcpp::NumericVector compute_Ddenom(Rcpp::NumericMatrix snpQ2, Rcpp::IntegerMatrix f2idx, Rcpp::LogicalVector verbose);
RcppExport SEXP _poolfstat_compute_Ddenom(SEXP snpQ2SEXP, SEXP f2idxSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type snpQ2(snpQ2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type f2idx(f2idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Ddenom(snpQ2, f2idx, verbose));
    return rcpp_result_gen;
END_RCPP
}
// compute_Q_bjmeans
Rcpp::NumericVector compute_Q_bjmeans(Rcpp::NumericMatrix snpQ, Rcpp::IntegerVector snp_bj_id, Rcpp::LogicalVector verbose);
RcppExport SEXP _poolfstat_compute_Q_bjmeans(SEXP snpQSEXP, SEXP snp_bj_idSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type snpQ(snpQSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type snp_bj_id(snp_bj_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Q_bjmeans(snpQ, snp_bj_id, verbose));
    return rcpp_result_gen;
END_RCPP
}
// compute_F2_bjmeans
Rcpp::NumericVector compute_F2_bjmeans(Rcpp::NumericMatrix snpQ1, Rcpp::NumericMatrix snpQ2, Rcpp::IntegerMatrix q1_idx, Rcpp::IntegerVector snp_bj_id, Rcpp::LogicalVector verbose);
RcppExport SEXP _poolfstat_compute_F2_bjmeans(SEXP snpQ1SEXP, SEXP snpQ2SEXP, SEXP q1_idxSEXP, SEXP snp_bj_idSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type snpQ1(snpQ1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type snpQ2(snpQ2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type q1_idx(q1_idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type snp_bj_id(snp_bj_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_F2_bjmeans(snpQ1, snpQ2, q1_idx, snp_bj_id, verbose));
    return rcpp_result_gen;
END_RCPP
}
// compute_Ddenom_bjmeans
Rcpp::NumericVector compute_Ddenom_bjmeans(Rcpp::NumericMatrix snpQ2, Rcpp::IntegerMatrix f2idx, Rcpp::IntegerVector snp_bj_id, Rcpp::LogicalVector verbose);
RcppExport SEXP _poolfstat_compute_Ddenom_bjmeans(SEXP snpQ2SEXP, SEXP f2idxSEXP, SEXP snp_bj_idSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type snpQ2(snpQ2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type f2idx(f2idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type snp_bj_id(snp_bj_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Ddenom_bjmeans(snpQ2, f2idx, snp_bj_id, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_poolfstat_scan_allele_info", (DL_FUNC) &_poolfstat_scan_allele_info, 1},
    {"_poolfstat_extract_vscan_counts", (DL_FUNC) &_poolfstat_extract_vscan_counts, 3},
    {"_poolfstat_extract_nonvscan_counts", (DL_FUNC) &_poolfstat_extract_nonvscan_counts, 4},
    {"_poolfstat_extract_allele_names", (DL_FUNC) &_poolfstat_extract_allele_names, 2},
    {"_poolfstat_compute_Ddenom", (DL_FUNC) &_poolfstat_compute_Ddenom, 3},
    {"_poolfstat_compute_Q_bjmeans", (DL_FUNC) &_poolfstat_compute_Q_bjmeans, 3},
    {"_poolfstat_compute_F2_bjmeans", (DL_FUNC) &_poolfstat_compute_F2_bjmeans, 5},
    {"_poolfstat_compute_Ddenom_bjmeans", (DL_FUNC) &_poolfstat_compute_Ddenom_bjmeans, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_poolfstat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
