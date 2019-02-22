#' Plot a region
#'
#' \code{plotRegion} provides the coverage of the replicates in a
#' \code{ggplot2} object.
#'
#' @param object An \code{srnadiff} object.
#' @param region A \code{GenomicRange} object.
#' @return A \code{ggplot2} object.
#' @examples
#' exp <- sRNADiffExample()
#' exp <- runAll(exp)
#' plotRegion(exp, regions(exp, 0.05)[1])
#'
#' @export
plotRegion <- function(object, region) {
    count     <- NULL
    pos       <- NULL
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
    sel       <- lapply(object@bamFiles, function (b) scanBam(b, param=param))
    cov       <- lapply(sel, function (c)
                    as.numeric(coverage(IRanges(c[["pos"]],
                        width=c[["qwidth"]]), width=e)[s:e]))
    nSamples  <- length(object@replicates)
    covTidy   <- data.frame(pos   =rep(seq(s, e), nSamples),
                            sample=rep(object@replicates, each=w),
                            count =c(apply(rbind(as.list(cov)), 1, unlist)))
    return(qplot(pos, count, data=covTidy, facets=sample~., geom='line',
            xlab=seqnames(region)[1], ylab='', main=names(region)[1]) +
            annotate("rect", xmin=sOr, xmax=eOr, ymin=-Inf, ymax=Inf, alpha=.2))
}
