#' Read the coverage and extract expressed regions
#'
#' @param object An \code{srnadiff} object.
#' @return       The selected values: a \code{list} of \code{vector}s of
#'                   \code{integer}s.
buildDataHmm <- function (object) {
    counts           <- rcpp_buildHmm(object@lengths, object@values,
                                        object@chromosomeSizes)
    counts           <- do.call(rbind, counts)
    colnames(counts) <- object@replicates
    rownames(counts) <- paste("data", seq(dim(counts)[1]), sep="_")
    return(counts)
}

#' Compute p-values of the selected counts.
#'
#' @param object An \code{srnadiff} object.
#' @param counts The counts: a \code{list} of \code{vector}s or \code{integer}s.
#' @return The p-values: a \code{list} of \code{numeric}.
computePvalues <- function(object, counts) {
    dds     <- DESeqDataSetFromMatrix(  countData=counts,
                                        colData  =object@design,
                                        design   =~ condition)
    dds     <- DESeq(dds)
    res     <- results(dds)
    pvalues <- res$padj
    rm(dds, res)
    gc()
    return(pvalues)
}

#' Initialize and run the HMM.
#'
#' @param object  An \code{srnadiff} object.
#' @param counts  The counts: a \code{list} of \code{vector}s or
#'                    \code{integer}s.
#' @param pvalues The p-values: a \code{list} of \code{numeric}.
#' @return A GRanges.
runHmm <- function(object, counts, pvalues) {
    NO_DIFF_THRESHOLD <- 0.85
    DIFF_THRESHOLD    <- 0.15
    NO_DIFF_TO_DIFF   <- 0.000001
    DIFF_TO_NO_DIFF   <- 0.001
    transitions       <- t(matrix(data=c(1 - NO_DIFF_TO_DIFF, NO_DIFF_TO_DIFF,
                            DIFF_TO_NO_DIFF, 1 - DIFF_TO_NO_DIFF), ncol=2))
    transitions       <- -log10(transitions / rowSums(transitions))
    emissions         <- t(matrix(data=c(1, 9, 9, 1), ncol=2))
    emissions         <- -log10(emissions / rowSums(emissions))
    starts            <- c(DIFF_TO_NO_DIFF / NO_DIFF_TO_DIFF,
                            NO_DIFF_TO_DIFF / DIFF_TO_NO_DIFF)
    starts            <- -log10(starts / sum(starts))
    ranges            <- rcpp_viterbi(object@chromosomeSizes, transitions,
                                        emissions, starts, counts, pvalues,
                                        object@lengths, object@values)
    return(GRanges(ranges))
}

#' Segmentation of the genome using an HMM.
#'
#' @param object An \code{srnadiff} object.
#' @return       A \code{GRanges} object.
runAllHmm <- function(object) {
    if (object@skipHmm) {
        return(GRanges())
    }
    message("Starting HMM step...")
    message("  Building data...")
    counts <- buildDataHmm(object)
    message("  ... data built")
    pvalues <- computePvalues(object, counts)
    message("  Running HMM...")
    intervals <- runHmm(object, counts, pvalues)
    message("  ... HMM run.")
    names(intervals) <- paste("hmm", seq(length(intervals)), sep="_")
    message(paste0(c("  ... ", length(intervals), " regions found.")))
    rm(counts, pvalues)
    gc()
    message("... HMM step done.")
    return(intervals)
}
