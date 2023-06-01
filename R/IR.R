##- Segmentation of the genome using the IR method ---------------------------#
##----------------------------------------------------------------------------#
runIR <- function(object) {

    message("Starting IR step...")

    logFC <- computeLogFoldChange(object)

    intervals <- rcpp_ir(logFC,
                        parameters(object)$minSize,
                        parameters(object)$maxSize,
                        parameters(object)$minLogFC)

    if (dim(intervals)[1] > 0) {
        intervals <- GRanges(intervals)
        names(intervals) <- paste(seqnames(intervals), start(intervals),
                                  end(intervals), sep="_")
    }

    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... IR step done.\n")
    return(intervals)
}
