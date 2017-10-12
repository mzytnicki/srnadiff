#' An S4 class to represent sRNA-Seq data for differential expression.
#'
#' @slot annotation      The annotation in GRanges format.
#' @slot bamFileNames    The name of one read file in BAM format.
#' @slot bamFiles        The BAM files in a \code{BamFileList}.
#' @slot chromosomes     The names of the chromosomes.
#' @slot chromosomeSizes The sizes of the chromosomes.
#' @slot replicates      The names of the replicates.
#' @slot conditions      The condition to which each replicate belongs.
#' @slot coverages       The coverages, a vector of \code{RLE}.
#' @slot lengths         The lengths parts of the coverages.
#' @slot values          The values parts of the coverages.
#' @slot design          Experimental design, a \code{DataFrame} for
#'                           \code{DESeq2}
#' @slot pValue          The adjusted p-value threshold
#' @slot skipAnnotation  Whether to skip the annotation strategy step
#' @slot skipNaive       Whether to skip the naive strategy step
#' @slot skipHmm         Whether to skip the HMM strategy step
setClass("sRNADiff",
            representation(
                annotation     ="GRanges",
                bamFileNames   ="vector",
                bamFiles       ="BamFileList",
                chromosomes    ="vector",
                chromosomeSizes="vector",
                replicates     ="vector",
                conditions     ="vector",
                coverages      ="vector",
                lengths        ="list",
                values         ="list",
                design         ="DataFrame",
                pValue         ="numeric",
                deseqData      ="DESeqDataSet",
                regions        ="GRanges",
                skipAnnotation ="logical",
                skipNaive      ="logical",
                skipHmm        ="logical",
                skipClustering ="logical"),
            prototype(
                bamFileNames=c(),
                replicates  =c(),
                conditions  =c()
            )
)

#' Constructor.
#'
#' @param annotation   The GRanges annotation
#' @param bamFileNames The name of one read file in BAM format.
#' @param replicates   The names of the replicates.
#' @param conditions   The condition to which each replicate belongs.
#' @param lazyload     Usual for S4 functions.
#' @return             An \code{sRNADiff} object.
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation  <- readWholeGenomeAnnotation(gtfFile)
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
#'
#' @export
sRNADiffExp <- function(annotation=NULL,
                        bamFileNames,
                        replicates,
                        conditions,
                        lazyload=FALSE) {
    message("Constructing object...")
    indexBam(bamFileNames)
    bamFiles <- BamFileList(bamFileNames, yieldSize=500000,
                            index=lapply(bamFileNames, paste0, ".bai"))
    names(bamFiles) <- tools::file_path_sans_ext(names(bamFiles))
    if (is.null(annotation)) {
        annotation <- GRanges()
    }
    object <- new("sRNADiff",
                    annotation     =annotation,
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
                    skipHmm        =FALSE,
                    skipClustering =FALSE
    )
    object@lengths <- lapply(object@chromosomes, function(chromosome)
        lapply(lapply(object@coverages, `[[`, chromosome), slot, "lengths"))
    object@values <- lapply(object@chromosomes, function(chromosome)
        lapply(lapply(object@coverages, `[[`, chromosome), slot, "values"))
    rownames(object@design) <- object@replicates
    message("... done.")
    return(object)
}
