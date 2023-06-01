##- Segmentation of the genome using the naive method ------------------------#
##----------------------------------------------------------------------------#
runNaive <- function(object) {

    message("Starting Naive step...")

    logFC <- abs(computeLogFoldChange(object))
    idReg <- (logFC >= parameters(object)$minLogFC)
    intervals <- lapply(idReg, IRanges)

    intervals <- lapply(intervals, function(IR, minSize, maxSize) {
                                        IR <- IR[(width(IR) >= minSize) &
                                                    (width(IR) <= maxSize)]
                                        return(IR) },
                        minSize=parameters(object)$minSize,
                        maxSize=parameters(object)$maxSize)

    intervals <- IRlist2GR(intervals)
    intervals <- reduce(intervals, min.gapwidth = parameters(object)$minGap)

    if (length(intervals) > 0) {
        names(intervals) <- paste(seqnames(intervals), start(intervals),
                                  end(intervals), sep="_")
    }

    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... Naive step done.\n")
    return(intervals)
}
