#' Plot a region
#'
#' @param object An \code{srnadiff} object.
#' @param region A \code{GenomicRange} object.
#' @return A \code{ggbio} object.
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
#' plotRegion(exp, diffRegions@regions[1])
#'
#' @export
plotRegion <- function(object, region) {
    cr        <- region
    sOr       <- start(region)[1]
    eOr       <- end(region)[1]
    wOr       <- eOr - sOr + 1
    s         <- max(1, sOr - 0.1 * wOr)
    e         <- eOr + 0.1 * wOr
    w         <- e - s + 1
    start(cr) <- s
    end(cr)   <- e
    param     <- ScanBamParam(which=cr, what=scanBamWhat())
    sel       <- sapply(object@bamFiles, function (b) scanBam(b, param=param))
    cov       <- lapply(sel, function (c)
        as.numeric(coverage(IRanges(c[["pos"]],
                                    width=c[["qwidth"]]), width=e)[s:e]))
    samples   <- sapply(strsplit(names(cov), "[.]"), "[[", 1)
    nSamples  <- length(samples)
    covTidy   <- data.frame(pos=rep(seq(s, e), nSamples),
                            sample=rep(samples, each=w),
                            count=c(apply(rbind(as.list(cov)), 1, unlist)))
    return(qplot(pos, count, data=covTidy, facets=sample~., geom='line',
                 xlab=seqnames(region)[1], ylab='',
                 main=names(region)[1]) +
               annotate("rect", xmin=sOr, xmax=eOr,
                        ymin=-Inf, ymax=Inf, alpha=.2))
}
