#' Initialize and run the slice method.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
runSlice <- function(object) {
    avgCounts <- lapply(split(object@coverages,
                                object@conditions == object@conditions[[1]]),
                            function (s) { round(Reduce('+', s) / length(s)) })
    lengths <- lapply(avgCounts, function (x) {lapply(x, slot, "lengths") })
    values  <- lapply(avgCounts, function (x) {lapply(x, slot, "values") })
    ranges  <- rcpp_slice(lengths, values, object@chromosomeSizes,
                            object@minDepth, object@minSize, object@maxSize,
                            object@minDifferences)
    if (length(ranges[[1]]) == 0) return(GRanges())
    return(GRanges(ranges))
}


#' Segmentation of the genome using a slice method.
#'
#' @param object An \code{srnadiff} object.
#' @return       A \code{GRanges} object.
runAllSlice <- function(object) {
    if (object@skipSlice) {
        return(GRanges())
    }
    message("Starting slice step...")
    intervals <- runSlice(object)
    if (length(intervals) > 0)
        names(intervals) <- paste("slice",seq(length(intervals)), sep="_")
    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... slice step done.")
    return(intervals)
}
