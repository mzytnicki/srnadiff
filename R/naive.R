build.model.naive <- function(bam.files) {
    mergeBam(bam.files, "tmp.bam", overwrite = TRUE)
    merged.reads <- readGAlignments("tmp.bam")
    merged.ranges <- granges(merged.reads)
    reduced.ranges <- reduce(merged.ranges, drop.empty.ranges = TRUE, ignore.strand = TRUE)
    names(reduced.ranges) <- as.character(reduced.ranges, ignore.strand = TRUE)
    return(reduced.ranges)
}

run.all.naive <- function(object) {
    annotation       <- build.model.naive(object@bam.files)
    counts           <- count.features(object, annotation)
    diff.data        <- run.deseq2(object, counts)
    return(rownames(diff.data))
}
