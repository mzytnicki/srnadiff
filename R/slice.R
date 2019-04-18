##- Segmentation of the genome using the slice method ------------------------#
##----------------------------------------------------------------------------#
runSlice <- function(object) {

    message("Starting Slice step...")

    logFC <- computeLogFoldChange(object)

    intervals <- rcpp_ir(logFC,
                        parameters(object)$minSize,
                        parameters(object)$maxSize,
                        parameters(object)$minLogFC)

    if (dim(intervals)[1] > 0) {
        intervals <- GRanges(intervals)
        names(intervals) <- paste("slice", seq(length(intervals)), sep="_")
    }

    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... Slice step done.\n")
    return(intervals)
}
