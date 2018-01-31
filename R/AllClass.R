#' An S4 class to represent sRNA-Seq data for differential expression.
#'
#' @slot annotation        The annotation in GRanges format.
#' @slot bamFileNames      The name of one read file in BAM format.
#' @slot bamFiles          The BAM files in a \code{BamFileList}.
#' @slot chromosomes       The names of the chromosomes.
#' @slot chromosomeSizes   The sizes of the chromosomes.
#' @slot replicates        The names of the replicates.
#' @slot conditions        The condition to which each replicate belongs.
#' @slot coverages         The coverages, a vector of \code{RLE}.
#' @slot lengths           The lengths parts of the coverages.
#' @slot values            The values parts of the coverages.
#' @slot design            Experimental design, a \code{DataFrame} for
#'                           \code{DESeq2}
#' @slot regions           A \code{GenomicRanges} of the possibly differentially
#'                           expressed region
#' @slot minDepth          Minimum depth to consider to find regions
#' @slot minSize           Minimum region size
#' @slot maxSize           Maximum region size
#' @slot mergeDistance     Distance to merge consecutive region
#' @slot minDifferences    Minimum number of different nt between two regions
#' @slot noDiffToDiff      Transition probability
#' @slot diffToNoDiff      Transition probability
#' @slot emission          Emission probability
#' @slot emissionThreshold Emission threshold
#' @slot skipAnnotation    Whether to skip the annotation strategy step
#' @slot skipNaive         Whether to skip the naive strategy step
#' @slot skipHmm           Whether to skip the HMM strategy step
#' @slot skipSlice         Whether to skip the slice strategy step
#' @slot nThreads          Number of threads
setClass("sRNADiff",
            representation(
                annotation       ="GRanges",
                bamFileNames     ="vector",
                bamFiles         ="BamFileList",
                chromosomes      ="vector",
                chromosomeSizes  ="vector",
                replicates       ="vector",
                conditions       ="vector",
                coverages        ="vector",
                lengths          ="list",
                values           ="list",
                design           ="DataFrame",
                regions          ="GRanges",
                minDepth         ="numeric",
                minSize          ="numeric",
                maxSize          ="numeric",
                mergeDistance    ="numeric",
                minDifferences   ="numeric",
                noDiffToDiff     ="numeric",
                diffToNoDiff     ="numeric",
                emission         ="numeric",
                emissionThreshold="numeric",
                minOverlap       ="numeric",
                skipAnnotation   ="logical",
                skipNaive        ="logical",
                skipHmm          ="logical",
                skipSlice        ="logical",
                nThreads         ="numeric"),
            prototype(
                annotation=GRanges()
            ),
            validity=function(object) {
                nBams  <- length(object@bamFileNames)
                nReps  <- length(object@replicates)
                nConds <- length(object@conditions)
                if (nBams != nReps) {
                    return(paste0("The number of input BAM files should be ",
                        "equal to the number of replicates (", nBams, " and ",
                        nReps, " resp.)."))
                }
                if (nBams != nConds) {
                    return(paste0("The number of input BAM files should be ",
                        "equal to the number of conditions (", nBams, " and ",
                        nConds, " resp.)."))
                }
                return(TRUE)
            }
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
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
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
    indexFileNames <- paste0(bamFileNames, ".bai")
    unIndexedFiles <- bamFileNames[! file.exists(indexFileNames)]
    if (length(unIndexedFiles) > 0) indexBam(unIndexedFiles)
    bamFiles <- BamFileList(lapply(bamFileNames,
                    function (b) {
                        BamFile(b, yieldSize=500000, paste0(b, ".bai"))}))
    names(bamFiles) <- tools::file_path_sans_ext(names(bamFiles))
    if (is.null(annotation)) {
        annotation <- GRanges()
    }
    object <- new("sRNADiff",
                    annotation       =annotation,
                    bamFileNames     =bamFileNames,
                    bamFiles         =bamFiles,
                    replicates       =replicates,
                    conditions       =conditions,
                    design           =DataFrame(condition=conditions),
                    chromosomes      =seqlevels(bamFiles[[1]]),
                    chromosomeSizes  =seqlengths(bamFiles[[1]]),
                    coverages        =lapply(bamFiles, coverage),
                    skipAnnotation   =FALSE,
                    skipNaive        =FALSE,
                    skipHmm          =FALSE,
                    skipSlice        =FALSE,
                    minDepth         =10,
                    minSize          =18,
                    maxSize          =1000000,
                    mergeDistance    =100,
                    minDifferences   =20,
                    noDiffToDiff     =0.001,
                    diffToNoDiff     =0.000001,
                    emission         =0.9,
                    emissionThreshold=0.1,
                    minOverlap       =10,
                    nThreads         =1
    )
    object@lengths <- lapply(object@chromosomes, function(chromosome)
        lapply(lapply(object@coverages, `[[`, chromosome), slot, "lengths"))
    object@values <- lapply(object@chromosomes, function(chromosome)
        lapply(lapply(object@coverages, `[[`, chromosome), slot, "values"))
    rownames(object@design) <- object@replicates
    message("... done.")
    return(object)
}


#' Example constructor
#' @return An \code{srnadiff} object
#'
#' @examples
#' exp <- sRNADiffExample()
#'
#' @export
sRNADiffExample <- function() {
    dir        <- system.file("extdata", package="srnadiff", mustWork=TRUE)
    data       <- read.csv(file.path(dir, "data.csv"))
    gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf.gz")
    annotation <- readWholeGenomeAnnotation(gtfFile)
    bamFiles   <- file.path(dir, data$FileName)
    replicates <- data$SampleName
    conditions <- factor(data$Condition)
    object     <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
    return(object)
}
