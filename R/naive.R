build.model.naive <- function(bam.files) {
    names(reduced.ranges) <- as.character(reduced.ranges, ignore.strand = TRUE)
}

run.all.naive <- function(object) {
    mergeBam(object@bam.files, "tmp.bam", overwrite = TRUE)
    merged.reads        <- readGAlignments("tmp.bam")
    merged.ranges       <- granges(merged.reads)
    reduced.ranges      <- reduce(merged.ranges, drop.empty.ranges = TRUE, ignore.strand = TRUE, with.revmap = TRUE)
    sizes               <- sapply(mcols(reduced.ranges)$revmap, length)
    sized.ranges        <- reduced.ranges[sizes >= 10]
    mcols(sized.ranges) <- NULL
    names(sized.ranges) <- paste("naive", seq(length(sized.ranges)), sep = "_")
    return(sized.ranges)
}
