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
setMethod("show",
          c(object = "sRNADiff"),
          definition = function (object) {
              cat(
                  "Object of class sRNADiff.\n",
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
                  "\np-value threshold: ",
                  object@pValue,
                  "\n"
                 , sep = "")
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
    allSets       <- GRangesList(c(setAnnotation, setNaive, setHmm))
    message("Computing differential expression...")
    counts <- do.call(rbind, lapply(list(setAnnotation, setNaive, setHmm),
                 function (set) { return(summarizeOverlaps(
                                           features     =set,
                                           reads        =object@bamFiles,
                                           mode         ="IntersectionNotEmpty",
                                           singleEnd    =TRUE,
                                           ignore.strand=TRUE,
                                           fragments    =FALSE))}))
    colnames(counts)           <- object@replicates
    colData(counts)            <- object@design
    dds                        <- DESeqDataSet(counts, design=~condition)
    dds                        <- dds[rowSums(counts(dds)) > 1, ]
    dds                        <- DESeq(dds)
    counts                     <- counts[rownames(counts) %in% rownames(dds)]
    regions                    <- counts@rowRanges
    mcols(regions)             <- as.data.frame(results(dds))
    padj                       <- mcols(regions)$padj
    regions                    <- regions[! is.na(padj)]
    dds                        <- dds[! is.na(padj)]
    padj                       <- padj[! is.na(padj)]
    sizes                      <- width(regions)
    overlaps                   <- findOverlaps(regions, regions)
    overlaps                   <- overlaps[overlaps@from != overlaps@to]
    dominance                  <- padj[overlaps@from] < padj[overlaps@to] |
        (padj[overlaps@from] == padj[overlaps@to] &
             sizes[overlaps@from] < sizes[overlaps@to]) |
        (padj[overlaps@from] == padj[overlaps@to] &
             sizes[overlaps@from] == sizes[overlaps@to] &
             overlaps@from < overlaps@to)
    regions                    <- regions[-overlaps@to[dominance]]
    names(regions)             <- paste0(seqnames(regions), "_", start(regions),
                                         "_", end(regions))
    deseqData                  <- dds[rownames(dds) %in% rownames(regions)]
    message("... done.")
    object@deseqData           <- dds
    object@regions             <- regions
    return(object)
}
