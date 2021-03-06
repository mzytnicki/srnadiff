# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_buildHmm <- function(coverages, minDepth) {
    .Call('_srnadiff_rcpp_buildHmm', PACKAGE = 'srnadiff', coverages, minDepth)
}

rcpp_viterbi <- function(coverages, transitions, emissions, emissionThreshold, starts, counts, pvalues, minDepth, minSize, maxSize) {
    .Call('_srnadiff_rcpp_viterbi', PACKAGE = 'srnadiff', coverages, transitions, emissions, emissionThreshold, starts, counts, pvalues, minDepth, minSize, maxSize)
}

rcpp_ir <- function(logFoldChanges, minLength, maxLength, minLFC) {
    .Call('_srnadiff_rcpp_ir', PACKAGE = 'srnadiff', logFoldChanges, minLength, maxLength, minLFC)
}

rcpp_normalization <- function(coverages, librarySizes) {
    .Call('_srnadiff_rcpp_normalization', PACKAGE = 'srnadiff', coverages, librarySizes)
}

