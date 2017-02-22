#' Segmentation of the genome in a naive way.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
runAllNaive <- function(object) {
    if (object@skipNaive) {
        return(GRanges())
    }
    tmpFileName        <- mergeBam(object@bamFiles,
                                   tempfile(pattern = "mergeTmp",
                                            fileext = ".bam"),
                                   overwrite=TRUE)
    mergedReads        <- readGAlignments(tmpFileName)
    mergedRanges       <- granges(mergedReads)
    reducedRanges      <- reduce(mergedRanges, drop.empty.ranges=TRUE,
                                  ignore.strand=TRUE, with.revmap=TRUE)
    sizes              <- sapply(mcols(reducedRanges)$revmap, length)
    sizedRanges        <- reducedRanges[sizes >= 10]
    mcols(sizedRanges) <- NULL
    names(sizedRanges) <- paste("naive", seq(length(sizedRanges)), sep="_")
    file.remove(tmpFileName)
    return(sizedRanges)
}
