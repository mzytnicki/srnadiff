#' Segmentation of the genome in a naive way.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
runAllNaive <- function(object) {
    if (object@skipNaive) {
        return(GRanges())
    }
    message("Starting Naive step...")
#    ranges <- rcpp_naive(object@lengths, object@values,
#                            object@chromosomeSizes, object@normalizationFactors,
#                            object@minDepth, object@mergeDistance,
#                            object@minSize)
#    if (length(ranges[[1]]) == 0) {
#        ranges <- GRanges()
#    }
#    else {
#        ranges <- GRanges(ranges)
#    }
    up <- as(object@logFC >= object@minLogFC, "GRanges")
    up <- up[up$score]
    up <- reduce(up, min.gapwidth=object@mergeDistance)
    down <- as(object@logFC <= -object@minLogFC, "GRanges")
    down <- down[down$score]
    down <- reduce(down, min.gapwidth=object@mergeDistance)
    ranges <- c(up, down)
    if (length(ranges) > 0) {
        ranges <- ranges[(width(ranges) >= object@minSize &
                              width(ranges) <= object@maxSize), ]
        mcols(ranges) <- NULL
        names(ranges) <- paste("naive", seq(length(ranges)), sep="_")
    }
    message(paste0(c("  ... ", length(ranges), " regions found.")))
    message("... Naive step done.")
    return(ranges)
}
