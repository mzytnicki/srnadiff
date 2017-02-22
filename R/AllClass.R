#' An S4 class to represent sRNA-Seq data for differential expression.
#'
#' @slot gtfFileName     The name of the annotation in GTF format.
#' @slot bamFileNames    The name of one read file in BAM format.
#' @slot bamFiles        The BAM files in a \code{BamFileList}.
#' @slot chromosomes     The names of the chromosomes.
#' @slot chromosomeSizes The sizes of the chromosomes.
#' @slot replicates      The names of the replicates.
#' @slot conditions      The condition to which each replicate belongs.
#' @slot coverages       The coverages, a vector of \code{RLE}.
#' @slot design          Experimental design, a \code{DataFrome} for \code{DESeq2}
#' @slot pValue          The adjusted p-value threshold
#' @slot skipAnnotation  Whether to skip the annotation strategy step
#' @slot skipNaive       Whether to skip the naive strategy step
#' @slot skipHmm         Whether to skip the HMM strategy step
setClass("sRNADiff",
         representation(
             gtfFileName    ="character",
             bamFileNames   ="vector",
             bamFiles       ="BamFileList",
             chromosomes    ="vector",
             chromosomeSizes="vector",
             replicates     ="vector",
             conditions     ="vector",
             coverages      ="vector",
             design         ="DataFrame",
             pValue         ="numeric",
             skipAnnotation ="logical",
             skipNaive      ="logical",
             skipHmm        ="logical"),
         prototype(
             gtfFileName =NA_character_,
             bamFileNames=c(),
             replicates  =c(),
             conditions  =c()))

#' Constructor.
#'
#' @param gtfFileName  The name of the annotation in GTF format.
#' @param bamFileNames The name of one read file in BAM format.
#' @param replicates   The names of the replicates.
#' @param conditions   The condition to which each replicate belongs.
#' @param lazyload     Usual for S4 functions.
#' @return             An \code{sRNADiff} object.
#' @export
sRNADiffExp <- function(gtfFileName=NA_character_,
                        bamFileNames,
                        replicates,
                        conditions,
                        lazyload=FALSE) {
    bamFiles <- BamFileList(bamFileNames)
    object <- new("sRNADiff",
                  gtfFileName    =gtfFileName,
                  bamFileNames   =bamFileNames,
                  bamFiles       =bamFiles,
                  replicates     =replicates,
                  conditions     =conditions,
                  design         =DataFrame(condition=conditions),
                  chromosomes    =seqlevels(bamFiles[[1]]),
                  chromosomeSizes=seqlengths(bamFiles[[1]]),
                  coverages      =sapply(bamFiles, coverage),
                  pValue         =0.05,
                  skipAnnotation =FALSE,
                  skipNaive      =FALSE,
                  skipHmm        =FALSE
    )
    rownames(object@design) <- object@replicates
    register(MulticoreParam(8))
    return(object)
}
