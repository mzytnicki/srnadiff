#' Segmentation of the genome using a IR method.
#'
#' @param object An \code{srnadiff} object.
#' @return       A \code{GRanges} object.
runAllIR <- function(object) {
    if (object@skipIR) {
        return(GRanges())
    }
    message("Starting IR step...")

    ranges <- rcpp_ir(object@logFC, object@minSize, object@maxSize,
                        object@minLogFC);
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
