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

    ranges <- rcpp_ir(object@logFC, object@minSize, object@maxSize, object@minLogFC);
    if (dim(ranges)[1] > 0) {
        ranges <- GRanges(ranges)
        names(ranges) <- paste("IR", seq(length(ranges)), sep="_")
        ranges$method <- as.factor("IR")
    }
    else {
        ranges <- GRanges()
    }

    message(paste0(c("  ... ", length(ranges), " regions found.")))
    message("... IR step done.")
    return(ranges)
}
