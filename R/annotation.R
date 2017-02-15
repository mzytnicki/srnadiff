run.all.annotation <- function(object) {
    database          <- makeTxDbFromGFF(object@gtf.file.name, format="gtf", circ_seqs=character())
    annotation        <- genes(database)
    mcols(annotation) <- NULL
    return (annotation)
}
