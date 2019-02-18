#' Initialize and run the IR method.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
#runSlice <- function(object) {
#    avgCounts <- lapply(split(object@coverages,
#                                object@conditions == object@conditions[[1]]),
#                            function (s) { round(Reduce('+', s) / length(s)) })
#    lengths <- lapply(avgCounts, function (x) {lapply(x, slot, "lengths") })
#    values  <- lapply(avgCounts, function (x) {lapply(x, slot, "values") })
#    ranges  <- rcpp_slice(lengths, values, object@chromosomeSizes,
#                            object@minDepth, object@minSize, object@maxSize,
#                            object@minDifferences)
#    if (length(ranges[[1]]) == 0) return(GRanges())
#    return(GRanges(ranges))
#}


#' Segmentation of the genome using a IR method.
#'
#' @param object An \code{srnadiff} object.
#' @return       A \code{GRanges} object.
runAllIR <- function(object) {
    if (object@skipIR) {
        return(GRanges())
    }
    message("Starting IR step...")

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
    intervals <- rcpp_ir(object@logFC, ranges, object@minSize, object@maxSize, object@minLogFC);
    intervals <- GRanges(intervals)
    if (length(intervals) > 0) {
        names(intervals) <- paste("IR", seq(length(intervals)), sep="_")
        intervals$method <- "IR"
    }

#    intervals <- runSlice(object)
#    if (length(intervals) > 0)
#        names(intervals) <- paste("slice",seq(length(intervals)), sep="_")
    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... IR step done.")
    return(intervals)
}
