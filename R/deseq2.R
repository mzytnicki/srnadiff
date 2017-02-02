count.features <- function(object, annotation) {
    counts   <- summarizeOverlaps(features      = annotation,
                                  reads         = object@bam.files,
                                  mode          = "IntersectionNotEmpty",
                                  singleEnd     = TRUE,
                                  ignore.strand = TRUE,
                                  fragments     = FALSE)
    colnames(counts) <- object@replicates
    colData(counts)  <- object@design
    return(counts)
}

run.deseq2 <- function(object, counts) {
    dds <- DESeqDataSet(counts, design = ~ condition)
    dds <- dds[ rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)
    res <- results(dds)
    res.05 <- res[ res$padj < 0.05 & ! is.na(res$padj), ]
    print(res)
    print(res.05)
    return(res.05)
}
