setClass("sRNADiff",
         representation(
             gtf.file.name        = "character",
             bam.file.names       = "vector",
             bam.files            = "BamFileList",
             chromosomes          = "vector",
             chromosome.sizes     = "vector",
             replicates           = "vector",
             conditions           = "vector",
             coverages            = "vector",
             design               = "DataFrame"),
         prototype(
             gtf.file.name        = NA_character_,
             bam.file.names       = c(),
             replicates           = c(),
             conditions           = c()))

sRNADiff <- function(gtf.file.name, bam.file.names, replicates, conditions, lazyload=FALSE) {
    bam.files <- BamFileList(bam.file.names)
    object <- new("sRNADiff",
                  gtf.file.name           = gtf.file.name,
                  bam.file.names          = bam.file.names,
                  bam.files               = bam.files,
                  replicates              = replicates,
                  conditions              = conditions,
                  design                  = DataFrame(condition = factor(conditions)),
                  chromosomes             = seqlevels(bam.files[[1]]),
                  chromosome.sizes        = seqlengths(bam.files[[1]]),
                  coverages               = sapply(bam.files, coverage)
    )
    rownames(object@design) = object@replicates
    register(MulticoreParam(8))
    return(object)
}
