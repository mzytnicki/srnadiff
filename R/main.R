#' Overloading the show method
#' @param object An \code{srnadiff} object.
#' @return       A description of the object.
#'
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation  <- readWholeGenomeAnnotation(gtfFile)
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
#' exp
#'
#' @export
setMethod(  f         ="show",
            signature ="sRNADiff",
            definition=function (object) {
                cat("Object of class sRNADiff.\n",
                    "Annotation (length = ",
                    length(object@annotation),
                    ")\nBAM files (length = ",
                    length(object@bamFileNames),
                    "):\n",
                    paste(' ', head(object@bamFileNames, 3)),
                    "...\nChromosomes (length = ",
                    length(object@chromosomes),
                    "):\n",
                    paste(' ', head(object@chromosomes, 3)),
                    "...\nReplicates (length = ",
                    length(object@replicates),
                    "):\n",
                    paste(' ', head(object@replicates, 3)),
                    "...\nConditions (length = ",
                    length(object@conditions),
                    "):\n",
                    paste(' ', head(object@conditions, 3)),
                    "...\nDifferentially expressed regions (length = ",
                    length(object@regions),
                    ")\nSkip annotation step: ",
                    object@skipAnnotation,
                    "\nSkip naive step: ",
                    object@skipNaive,
                    "\nSkip HMM step: ",
                    object@skipHmm,
                    "\nSkip clustering step: ",
                    object@skipClustering,
                    "\np-value threshold: ",
                    object@pValue,
                    "\n", sep = "")
            }
)

#' Get the output regions
#' @rdname regions-method
#' @param  object  An \code{srnadiff} object.
#' @return         The selected regions
#'
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation  <- readWholeGenomeAnnotation(gtfFile)
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
#' regions(exp)
#'
#' @export
setGeneric( name="regions",
            def =function(object) {
                standardGeneric("regions")
            }
)

#' @rdname  regions-method
#' @export
setMethod(  f        ="regions",
            signature= "sRNADiff",
            definition=function(object) {
                return(object@regions)
            }
)

#' Run the segmentation using 3 different methods, and reconciliate them.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation  <- readWholeGenomeAnnotation(gtfFile)
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(annotation, bamFiles, replicates, conditions)
#' diffRegions <- runAll(exp)
#'
#' @export
runAll <- function(object) {
    setAnnotation <- runAllAnnotation(object)
    setNaive      <- runAllNaive(object)
    setHmm        <- runAllHmm(object)
    setClustering <- runAllClustering(object)
    allSets       <- do.call("c", list(setAnnotation, setNaive, setHmm,
                        setClustering))
    message("Computing differential expression...")
    allSetsDT           <- as.data.frame(allSets)
    names(allSetsDT)    <- c("Chr", "Start", "End", "Width", "Strand")
    allSetsDT$GeneId    <- names(allSets)
    counts <- featureCounts(object@bamFileNames,
                            annot.ext=allSetsDT,
                            allowMultiOverlap=TRUE,
                            countMultiMappingReads=TRUE)
    counts              <- as.data.frame(counts$counts)
    rownames(counts)    <- names(allSets)
    colnames(counts)    <- object@replicates
    dds                 <- DESeqDataSetFromMatrix(  countData=counts,
                                                    colData  =object@design,
                                                    design   =~condition)
    names(dds)          <- names(allSets)
    dds                 <- dds[rowSums(counts(dds)) > 1, ]
    dds                 <- DESeq(dds)
    regions             <- allSets[names(allSets) %in% names(dds),]
    mcols(regions)      <- results(dds)
    padj                <- mcols(regions)$padj
    regions             <- regions[! is.na(padj)]
    dds                 <- dds[! is.na(padj)]
    padj                <- padj[! is.na(padj)]
    sizes               <- width(regions)
    overlaps            <- findOverlaps(regions, regions)
    overlaps            <- overlaps[overlaps@from != overlaps@to]
    if (length(overlaps) > 0) {
        dominance <- padj[overlaps@from] < padj[overlaps@to] |
            (padj[overlaps@from] == padj[overlaps@to] &
                sizes[overlaps@from] < sizes[overlaps@to]) |
            (padj[overlaps@from] == padj[overlaps@to] &
                sizes[overlaps@from] == sizes[overlaps@to] &
                overlaps@from < overlaps@to)
        regions   <- regions[-overlaps@to[dominance]]
    }
    names(regions)      <- paste0(seqnames(regions), "_", start(regions),
                                                            "_", end(regions))
    deseqData           <- dds[rownames(dds) %in% rownames(regions)]
    message("... done.")
    object@deseqData    <- dds
    object@regions      <- regions
    return(object)
}
