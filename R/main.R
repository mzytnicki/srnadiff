#' Run the segmentation using 3 different methods, and reconciliate them.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
#' @examples
#' dir         <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' data        <- read.csv(file.path(dir, "data.csv"))
#' gtfFile     <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' bamFiles    <- file.path(dir, data$FileName)
#' replicates  <- data$SampleName
#' conditions  <- factor(data$Condition)
#' exp         <- sRNADiffExp(gtfFile, bamFiles, replicates, conditions)
#' diffRegions <- runAll(exp)
#'
#' @export
runAll <- function(object) {
    setAnnotation              <- runAllAnnotation(object)
    setNaive                   <- runAllNaive(object)
    setHmm                     <- runAllHmm(object)
    allSets                    <- c(setAnnotation, setNaive, setHmm)
    allCounts                  <- countFeatures(object, allSets)
    allPvalues                 <- runDeseq2(object, allCounts)
    diffSets                   <- allSets[names(allSets) %in%
                                              rownames(allPvalues)]
    padj                       <- allPvalues$padj
    overlaps                   <- findOverlaps(diffSets, diffSets)
    undominated                <- rep(TRUE, length(overlaps))
    undominated[overlaps@from] <- (padj[overlaps@from] <= padj[overlaps@to])
    undominated[overlaps@to]   <- undominated[overlaps@to] ||
                                    (padj[overlaps@from] <= padj[overlaps@to])
    diffSets                   <- diffSets[undominated]
    print(diffSets)
    return(diffSets)
}
