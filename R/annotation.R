#' Parse an annotation file
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
runAllAnnotation <- function(object) {
    if (is.na(object@gtfFileName) || object@skipAnnotation) {
        return(GRanges())
    }
    database          <- makeTxDbFromGFF(object@gtfFileName,
                                         format="gtf",
                                         circ_seqs=character())
    annotation        <- genes(database)
    mcols(annotation) <- NULL
    return (annotation)
}
