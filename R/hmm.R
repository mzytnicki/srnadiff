#' Read the coverage and extract expressed regions
#'
#' @param object An \code{srnadiff} object.
#' @return       The selected values: a \code{list} of \code{vector}s of
#'                   \code{integer}s.
buildDataHmm <- function (object) {
    counts           <- rcpp_buildHmm(object@coverages, object@minDepth)
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
    dds     <- suppressMessages(DESeq(dds, parallel=(object@nThreads > 1)))
    res     <- results(dds)
    pvalues <- res$padj
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
    transitions       <- t(matrix(data=c(1 - object@noDiffToDiff,
                            object@noDiffToDiff, object@diffToNoDiff,
                            1 - object@diffToNoDiff), ncol=2))
    transitions       <- -log10(transitions / rowSums(transitions))
    emissions         <- t(matrix(data=c(1 - object@emission, object@emission,
                                            object@emission,
                                            1 - object@emission), ncol=2))
    emissions         <- -log10(emissions / rowSums(emissions))
    starts            <- c(object@noDiffToDiff / object@diffToNoDiff,
                            object@diffToNoDiff / object@noDiffToDiff)
    starts            <- -log10(starts / sum(starts))
    ranges            <- rcpp_viterbi(object@coverages, transitions,
                                        emissions, object@emissionThreshold,
                                        starts, counts, pvalues,
                                        object@minDepth, object@minSize,
                                        object@maxSize)
    if (nrow(ranges) == 0) return(GRanges())
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
    message("  Computing p-values...")
    pvalues <- computePvalues(object, counts)
    message("  ... values computed")
    message("  Running HMM...")
    intervals <- runHmm(object, counts, pvalues)
    message("  ... HMM run.")
    if (length(intervals) > 0) {
        names(intervals) <- paste("hmm", seq(length(intervals)), sep="_")
        intervals$method <- as.factor("HMM")
    }
    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... HMM step done.")
    return(intervals)
}
