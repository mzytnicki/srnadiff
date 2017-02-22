#' Count the number of reads in each condition in the given annotation.
#'
#' @param object An \code{srnadiff} object.
#' @param annotation A \code{GRanges} object.
#' @return A \code{list} of \code{vector} of \code{integer}s.
countFeatures <- function(object, annotation) {
    counts   <- summarizeOverlaps(features     =annotation,
                                  reads        =object@bamFiles,
                                  mode         ="IntersectionNotEmpty",
                                  singleEnd    =TRUE,
                                  ignore.strand=TRUE,
                                  fragments    =FALSE)
    colnames(counts) <- object@replicates
    colData(counts)  <- object@design
    return(counts)
}

#' Find differentially expressed features.
#'
#' @param object An \code{srnadiff} object.
#' @param counts A \code{list} of \code{vector} of \code{integer}s.
#' @return A \code{list} of feature names.
runDeseq2 <- function(object, counts) {
    dds   <- DESeqDataSet(counts, design=~condition)
    dds   <- dds[ rowSums(counts(dds)) > 1, ]
    dds   <- DESeq(dds)
    res   <- results(dds)
    res05 <- res[ res$padj < object@pValue & ! is.na(res$padj), ]
    return(res05)
}
