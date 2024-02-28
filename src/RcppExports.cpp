// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// barcodeMatch
DataFrame barcodeMatch(std::vector<std::string> seq, std::vector<std::string> barcodes, double mu, double sigma, double sigma_start, int k, int batch, int top, double cos_thresh, double alpha, int edit_thresh, const int UMI_len, const int flank);
RcppExport SEXP _Longcellsrc_barcodeMatch(SEXP seqSEXP, SEXP barcodesSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP sigma_startSEXP, SEXP kSEXP, SEXP batchSEXP, SEXP topSEXP, SEXP cos_threshSEXP, SEXP alphaSEXP, SEXP edit_threshSEXP, SEXP UMI_lenSEXP, SEXP flankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type barcodes(barcodesSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_start(sigma_startSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< int >::type top(topSEXP);
    Rcpp::traits::input_parameter< double >::type cos_thresh(cos_threshSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type edit_thresh(edit_threshSEXP);
    Rcpp::traits::input_parameter< const int >::type UMI_len(UMI_lenSEXP);
    Rcpp::traits::input_parameter< const int >::type flank(flankSEXP);
    rcpp_result_gen = Rcpp::wrap(barcodeMatch(seq, barcodes, mu, sigma, sigma_start, k, batch, top, cos_thresh, alpha, edit_thresh, UMI_len, flank));
    return rcpp_result_gen;
END_RCPP
}
// NeighborExtract
std::vector<std::vector<std::string> > NeighborExtract(std::vector<std::string> reads, std::vector<int> id, std::vector<int> pos, const int UMI_len, const int flank, const int bar_len);
RcppExport SEXP _Longcellsrc_NeighborExtract(SEXP readsSEXP, SEXP idSEXP, SEXP posSEXP, SEXP UMI_lenSEXP, SEXP flankSEXP, SEXP bar_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type reads(readsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type id(idSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type pos(posSEXP);
    Rcpp::traits::input_parameter< const int >::type UMI_len(UMI_lenSEXP);
    Rcpp::traits::input_parameter< const int >::type flank(flankSEXP);
    Rcpp::traits::input_parameter< const int >::type bar_len(bar_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(NeighborExtract(reads, id, pos, UMI_len, flank, bar_len));
    return rcpp_result_gen;
END_RCPP
}
// isos2exonids_index
List isos2exonids_index(std::vector<std::string> isoform, std::vector<int> start, std::vector<int> end, std::vector<std::string> exon_id, int mid_bias, int end_bias, int end_overlap, std::string nonsense_label, std::string split, std::string sep);
RcppExport SEXP _Longcellsrc_isos2exonids_index(SEXP isoformSEXP, SEXP startSEXP, SEXP endSEXP, SEXP exon_idSEXP, SEXP mid_biasSEXP, SEXP end_biasSEXP, SEXP end_overlapSEXP, SEXP nonsense_labelSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type isoform(isoformSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type end(endSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type exon_id(exon_idSEXP);
    Rcpp::traits::input_parameter< int >::type mid_bias(mid_biasSEXP);
    Rcpp::traits::input_parameter< int >::type end_bias(end_biasSEXP);
    Rcpp::traits::input_parameter< int >::type end_overlap(end_overlapSEXP);
    Rcpp::traits::input_parameter< std::string >::type nonsense_label(nonsense_labelSEXP);
    Rcpp::traits::input_parameter< std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(isos2exonids_index(isoform, start, end, exon_id, mid_bias, end_bias, end_overlap, nonsense_label, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// extractReads
DataFrame extractReads(std::vector<std::string> seq, std::vector<std::string> cigar, NumericVector pos, NumericMatrix annotation, const std::string strand, int toolkit, int end_flank, int splice_site_bin);
RcppExport SEXP _Longcellsrc_extractReads(SEXP seqSEXP, SEXP cigarSEXP, SEXP posSEXP, SEXP annotationSEXP, SEXP strandSEXP, SEXP toolkitSEXP, SEXP end_flankSEXP, SEXP splice_site_binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type cigar(cigarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type annotation(annotationSEXP);
    Rcpp::traits::input_parameter< const std::string >::type strand(strandSEXP);
    Rcpp::traits::input_parameter< int >::type toolkit(toolkitSEXP);
    Rcpp::traits::input_parameter< int >::type end_flank(end_flankSEXP);
    Rcpp::traits::input_parameter< int >::type splice_site_bin(splice_site_binSEXP);
    rcpp_result_gen = Rcpp::wrap(extractReads(seq, cigar, pos, annotation, strand, toolkit, end_flank, splice_site_bin));
    return rcpp_result_gen;
END_RCPP
}
// cigarProcess
List cigarProcess(std::string cigar);
RcppExport SEXP _Longcellsrc_cigarProcess(SEXP cigarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type cigar(cigarSEXP);
    rcpp_result_gen = Rcpp::wrap(cigarProcess(cigar));
    return rcpp_result_gen;
END_RCPP
}
// seqEnd
int seqEnd(int start, StringVector mark, NumericVector count);
RcppExport SEXP _Longcellsrc_seqEnd(SEXP startSEXP, SEXP markSEXP, SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< StringVector >::type mark(markSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type count(countSEXP);
    rcpp_result_gen = Rcpp::wrap(seqEnd(start, mark, count));
    return rcpp_result_gen;
END_RCPP
}
// isos_dis
DataFrame isos_dis(const std::vector<std::string> isoforms, const int thresh, const std::string split, const std::string sep);
RcppExport SEXP _Longcellsrc_isos_dis(SEXP isoformsSEXP, SEXP threshSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type isoforms(isoformsSEXP);
    Rcpp::traits::input_parameter< const int >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< const std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(isos_dis(isoforms, thresh, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// size_filter_cpp
NumericVector size_filter_cpp(NumericVector size, double ratio);
RcppExport SEXP _Longcellsrc_size_filter_cpp(SEXP sizeSEXP, SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(size_filter_cpp(size, ratio));
    return rcpp_result_gen;
END_RCPP
}
// splice_site_table_cpp
List splice_site_table_cpp(std::vector<std::string> isoform, const std::string split, const std::string sep, const int splice_site_thresh);
RcppExport SEXP _Longcellsrc_splice_site_table_cpp(SEXP isoformSEXP, SEXP splitSEXP, SEXP sepSEXP, SEXP splice_site_threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type isoform(isoformSEXP);
    Rcpp::traits::input_parameter< const std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< const std::string >::type sep(sepSEXP);
    Rcpp::traits::input_parameter< const int >::type splice_site_thresh(splice_site_threshSEXP);
    rcpp_result_gen = Rcpp::wrap(splice_site_table_cpp(isoform, split, sep, splice_site_thresh));
    return rcpp_result_gen;
END_RCPP
}
// extractTagFastq
DataFrame extractTagFastq(const char* fastq_path, const char* out_path, const std::string adapter, const int toolkit, const int window, const int step, const int len, const int polyA_bin, const int polyA_base_count, const int polyA_len);
RcppExport SEXP _Longcellsrc_extractTagFastq(SEXP fastq_pathSEXP, SEXP out_pathSEXP, SEXP adapterSEXP, SEXP toolkitSEXP, SEXP windowSEXP, SEXP stepSEXP, SEXP lenSEXP, SEXP polyA_binSEXP, SEXP polyA_base_countSEXP, SEXP polyA_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type fastq_path(fastq_pathSEXP);
    Rcpp::traits::input_parameter< const char* >::type out_path(out_pathSEXP);
    Rcpp::traits::input_parameter< const std::string >::type adapter(adapterSEXP);
    Rcpp::traits::input_parameter< const int >::type toolkit(toolkitSEXP);
    Rcpp::traits::input_parameter< const int >::type window(windowSEXP);
    Rcpp::traits::input_parameter< const int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const int >::type polyA_bin(polyA_binSEXP);
    Rcpp::traits::input_parameter< const int >::type polyA_base_count(polyA_base_countSEXP);
    Rcpp::traits::input_parameter< const int >::type polyA_len(polyA_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(extractTagFastq(fastq_path, out_path, adapter, toolkit, window, step, len, polyA_bin, polyA_base_count, polyA_len));
    return rcpp_result_gen;
END_RCPP
}
// baseCount
int baseCount(std::string seq, char base);
RcppExport SEXP _Longcellsrc_baseCount(SEXP seqSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< char >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(baseCount(seq, base));
    return rcpp_result_gen;
END_RCPP
}
// shareNeighbor
DataFrame shareNeighbor(std::vector<std::string> index, List neighbor, NumericVector count);
RcppExport SEXP _Longcellsrc_shareNeighbor(SEXP indexSEXP, SEXP neighborSEXP, SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< List >::type neighbor(neighborSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type count(countSEXP);
    rcpp_result_gen = Rcpp::wrap(shareNeighbor(index, neighbor, count));
    return rcpp_result_gen;
END_RCPP
}
// isoform2sites
std::vector<int> isoform2sites(std::string iso, const std::string split, const std::string sep);
RcppExport SEXP _Longcellsrc_isoform2sites(SEXP isoSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type iso(isoSEXP);
    Rcpp::traits::input_parameter< const std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< const std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(isoform2sites(iso, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// sites_chop
std::vector<int> sites_chop(std::vector<int> sites, int start, int end);
RcppExport SEXP _Longcellsrc_sites_chop(SEXP sitesSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type sites(sitesSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(sites_chop(sites, start, end));
    return rcpp_result_gen;
END_RCPP
}
// iso2_mid_dist
int iso2_mid_dist(std::string a, std::string b, const std::string split, const std::string sep);
RcppExport SEXP _Longcellsrc_iso2_mid_dist(SEXP aSEXP, SEXP bSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< const std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< const std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(iso2_mid_dist(a, b, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// iso2_mid_diff
NumericVector iso2_mid_diff(std::string a, std::string b, const int end_bias, const std::string split, const std::string sep);
RcppExport SEXP _Longcellsrc_iso2_mid_diff(SEXP aSEXP, SEXP bSEXP, SEXP end_biasSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type end_bias(end_biasSEXP);
    Rcpp::traits::input_parameter< const std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< const std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(iso2_mid_diff(a, b, end_bias, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// isoset_mid_diff
NumericMatrix isoset_mid_diff(std::vector<std::string> iso_set1, std::vector<std::string> iso_set2, const int thresh, const double overlap_thresh, const int end_bias, std::string split, std::string sep);
RcppExport SEXP _Longcellsrc_isoset_mid_diff(SEXP iso_set1SEXP, SEXP iso_set2SEXP, SEXP threshSEXP, SEXP overlap_threshSEXP, SEXP end_biasSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type iso_set1(iso_set1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type iso_set2(iso_set2SEXP);
    Rcpp::traits::input_parameter< const int >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const double >::type overlap_thresh(overlap_threshSEXP);
    Rcpp::traits::input_parameter< const int >::type end_bias(end_biasSEXP);
    Rcpp::traits::input_parameter< std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(isoset_mid_diff(iso_set1, iso_set2, thresh, overlap_thresh, end_bias, split, sep));
    return rcpp_result_gen;
END_RCPP
}
// umi_graph_table
std::vector<std::vector<int> > umi_graph_table(std::vector<std::string> umi, std::vector<std::string> isoform, std::vector<int> count, int sim_thresh, int iso_thresh, std::string split, std::string sep);
RcppExport SEXP _Longcellsrc_umi_graph_table(SEXP umiSEXP, SEXP isoformSEXP, SEXP countSEXP, SEXP sim_threshSEXP, SEXP iso_threshSEXP, SEXP splitSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type umi(umiSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type isoform(isoformSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type count(countSEXP);
    Rcpp::traits::input_parameter< int >::type sim_thresh(sim_threshSEXP);
    Rcpp::traits::input_parameter< int >::type iso_thresh(iso_threshSEXP);
    Rcpp::traits::input_parameter< std::string >::type split(splitSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(umi_graph_table(umi, isoform, count, sim_thresh, iso_thresh, split, sep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Longcellsrc_barcodeMatch", (DL_FUNC) &_Longcellsrc_barcodeMatch, 13},
    {"_Longcellsrc_NeighborExtract", (DL_FUNC) &_Longcellsrc_NeighborExtract, 6},
    {"_Longcellsrc_isos2exonids_index", (DL_FUNC) &_Longcellsrc_isos2exonids_index, 10},
    {"_Longcellsrc_extractReads", (DL_FUNC) &_Longcellsrc_extractReads, 8},
    {"_Longcellsrc_cigarProcess", (DL_FUNC) &_Longcellsrc_cigarProcess, 1},
    {"_Longcellsrc_seqEnd", (DL_FUNC) &_Longcellsrc_seqEnd, 3},
    {"_Longcellsrc_isos_dis", (DL_FUNC) &_Longcellsrc_isos_dis, 4},
    {"_Longcellsrc_size_filter_cpp", (DL_FUNC) &_Longcellsrc_size_filter_cpp, 2},
    {"_Longcellsrc_splice_site_table_cpp", (DL_FUNC) &_Longcellsrc_splice_site_table_cpp, 4},
    {"_Longcellsrc_extractTagFastq", (DL_FUNC) &_Longcellsrc_extractTagFastq, 10},
    {"_Longcellsrc_baseCount", (DL_FUNC) &_Longcellsrc_baseCount, 2},
    {"_Longcellsrc_shareNeighbor", (DL_FUNC) &_Longcellsrc_shareNeighbor, 3},
    {"_Longcellsrc_isoform2sites", (DL_FUNC) &_Longcellsrc_isoform2sites, 3},
    {"_Longcellsrc_sites_chop", (DL_FUNC) &_Longcellsrc_sites_chop, 3},
    {"_Longcellsrc_iso2_mid_dist", (DL_FUNC) &_Longcellsrc_iso2_mid_dist, 4},
    {"_Longcellsrc_iso2_mid_diff", (DL_FUNC) &_Longcellsrc_iso2_mid_diff, 5},
    {"_Longcellsrc_isoset_mid_diff", (DL_FUNC) &_Longcellsrc_isoset_mid_diff, 7},
    {"_Longcellsrc_umi_graph_table", (DL_FUNC) &_Longcellsrc_umi_graph_table, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_Longcellsrc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
