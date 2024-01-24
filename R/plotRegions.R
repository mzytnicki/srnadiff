#' Plot the coverage information surrounding genomic regions
#'
#' This function plot the coverage information surrounding genomic regions
#' while summarizing the annotation.
#'
#' This function provides a flexible genomic visualization framework by
#' displaying tracks in the sense of the \code{Gviz} package. Given
#' a region (or regions), four separate tracks are represented:
#' (1) \code{\link[Gviz]{GenomeAxisTrack}}, a horizontal axis with genomic
#' coordinate tickmarks for reference location to the displayed genomic
#' regions;
#' (2) \code{\link[Gviz]{GeneRegionTrack}}, if the \code{annot} argument is
#' passed, a track displaying all gene and/or sRNA annotation information
#' in a particular region;
#' (3) \code{\link[Gviz]{AnnotationTrack}}, regions are plotted as simple
#' boxes if no strand information is available, or as arrows to indicate
#' their direction; and
#' (4) \code{\link[Gviz]{DataTrack}}, plot the sample coverages surrounding
#' the genomic regions.
#'
#' The sample coverages can be plotted in various different forms as well
#' as combinations thereof. Supported plotting types are:
#'
#' \describe{
#' \item{}{\code{p:} simple dot plot.}
#' \item{}{\code{l:} lines plot.}
#' \item{}{\code{b:} combination of dot and lines plot.}
#' \item{}{\code{a:} lines plot of the sample-groups average (i.e., mean)
#'         values.}
#' \item{}{\code{confint:} confidence intervals for average values. In
#'         combination with \code{a} type.}
#' }
#'
#' @param object       An object of class \code{\link{srnadiffExp}} as returned
#'                     by \code{\link{srnadiff}} function.
#' @param region       A \code{GRanges} object. By example, ranges in the
#'                     output from \code{\link{regions}} used on an object of
#'                     class \code{\link{srnadiffExp}}.
#' @param normCvg      Boolean. If \code{TRUE} (the default), normalized
#'                     coverages are displayed, else, the raw coverages are
#'                     used.
#' @param annot        A \code{GRanges} or \code{data.frame} object containing
#'                     gene and/or sRNA annotation information.
#' @param allSignReg   Boolean. If \code{TRUE}, all differeltially expressed
#'                     regions contained in the annotated regions will be
#'                     displayed.
#' @param featureAnnot Character scalar. Feature annotation to be used, it
#'                     will be one from column names of \code{mcols(annot)}
#'                     if \code{annot} is a \code{GRanges} or one from
#'                     \code{colnames(annot)} if \code{annot} is a
#'                     \code{data.frame}.
#' @param flankReg     Integer value. If \code{flankReg} is positive, the
#'                     displayed genomic interval is extended by
#'                     \code{flankReg} nucleotides from the minimum
#'                     \code{start} and maximum \code{end} from the displayed
#'                     (ranges) regions.
#' @param colGroup     Character vector of length 2. The fill colors to be used
#'                     to indicate samples per group.
#' @param fillReg      Character vector of length 2. The fill colors to be used
#'                     to indicate 'up' or 'down' regulated regions.
#' @param fillAnnot    Character scalar. The fill color to be used for
#'                     annotated regions.
#' @param type         Character vector. The plot type, one of
#'                     (\code{"a", "l", "p", "b", "confint"}) or combinations
#'                     thereof. See Details for more information on the
#'                     individual plotting types.
#' @param trNames      Title for the tracks. By default \code{"DER"} and
#'                     \code{"coverage"}.
#' @param chrTitle     Boolean or character vector of length one. Defaults
#'                     \code{TRUE}, the chromosome name is added to the plot.
#'                     If character, this will be used for title.
#' @param legend       Boolean triggering the addition of a legend to the
#'                     (coverage) data track to indicate groups.
#' @param ...          Additional display parameters to control the look and
#'                     feel of the plots. See the "Display Parameters" section
#'                     for functions \code{\link{GenomeAxisTrack}},
#'                     \code{\link{GeneRegionTrack}},
#'                     \code{\link{AnnotationTrack}}, \code{\link{DataTrack}}
#'                     in the \code{Gviz} package.
#' @return             A list of GenomeGraph track objects, all inheriting
#'                     from class \code{GdObject}.
#'
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#'
#' plotRegions(srnaExp, regions(srnaExp)[1])
#' plotRegions(srnaExp, regions(srnaExp)[1], type = c("a", "confint"))
#'
#' @export
plotRegions <- function(object, region,
                        normCvg = TRUE,
                        annot = NULL,
                        allSignReg = TRUE,
                        featureAnnot = NULL,
                        flankReg = 10,
                        colGroup = c("#0080ff", "#ff00ff"),
                        fillReg = c("darkgreen", "darkred"),
                        fillAnnot = "#FFD58A",
                        type = "a",
                        chrTitle = TRUE,
                        trNames = c("DER", "coverage"),
                        legend = TRUE,
                        ...) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- internal function for character color checking -----------#
    ##------------------------------------------------------------#
    isColor <- function(x) { vapply(x, function(x) {
        tryCatch(is.matrix(col2rgb(x)), error=function(e) FALSE) },
        TRUE)
    }
    ##------------------------------------------------------------#

    ##- object
    if (!is(object, "srnadiffExp")) {
        stop("'object' must be an object of class 'srnadiffExp'.", call.=FALSE)
    }

    ##- region
    if (!is(region, "GRanges")) {
        stop("'region' must be an object of class 'GRanges'.", call.=FALSE)
    }

    ##- normCvg
    if (length(normCvg) != 1 || !is.logical(normCvg)) {
        stop("'normCvg' must be either TRUE or FALSE.", call.=FALSE)
    }

    ##- annot
    if (!is.null(annot)) {
        if (!all(class(annot) %in% c("GRanges", "data.frame"))) {
            stop("'annot' must be an object of class either 'GRanges'",
                 " or 'data.frame.", call.=FALSE)
        }

        if (is(annot, "data.frame")) {
            annot <- makeGRangesFromDataFrame(annot,
                                                keep.extra.columns = TRUE)
        }
    }

    ##- allSignReg
    if (length(allSignReg) != 1 | !is.logical(allSignReg)) {
        stop("'allSignReg' must be either TRUE or FALSE.", call.=FALSE)
    }

    ##- featureAnnot
    if (!is.null(annot)) {
        if (is.null(featureAnnot)) {
            featureAnnot <- colnames(mcols(annot))
            if (length(featureAnnot) > 0) {
                stop("'featureAnnot' should be one of {",
                     paste0(featureAnnot, collapse = ", "),
                     "}.", call.=FALSE)
            } else {
                warning("Non feature annotation in 'annot'.",
                        " 'annot' will be ignored.", call.=FALSE)
            }
        }
    }

    ##- flankReg
    if (length(flankReg) != 1) {
        stop("'flankReg' must be of length 1.", call.=FALSE)
    }

    if (is.null(flankReg) || !is.numeric(flankReg) || (flankReg < 0) ||
        !is.finite(flankReg)) {
        stop("invalid value for 'flankReg'. It must be a not negative",
             " integer.", call.=FALSE)
    }

    dec <- flankReg - trunc(flankReg)

    if (dec > 0) {
        stop("invalid value for 'flankReg'. It must be a not negative",
             " integer.", call.=FALSE)
    }

    ##- colGroup
    if (length(colGroup) != 2) {
        stop("'colGroup' must be a character vector of color names of",
            " length 2.", call.=FALSE)
    } else {
        if (any(!isColor(colGroup))) {
            stop("'colGroup' must be a character vector of recognized",
                 " colors.", call.=FALSE)
        }
    }

    ##- fillReg
    if (length(fillReg) != 2) {
        stop("'fillReg' must be a character vector of color names of",
            " length 2.", call.=FALSE)
    } else {
        if (any(!isColor(fillReg))) {
            stop("'fillReg' must be a character vector of recognized",
                 " colors.", call.=FALSE)
        }
    }

    ##- fillAnnot
    if (length(fillAnnot) != 1) {
        stop("'fillAnnot' must be a character vector of color names of",
            " length 1.", call.=FALSE)
    } else {
        if (any(!isColor(fillAnnot))) {
            stop("'fillAnnot' must be a character of recognized",
                 " colors.", call.=FALSE)
        }
    }

    ##- type
    choices <- c("a", "l", "p", "b", "confint")
    type <- choices[pmatch(type, choices)]

    if (any(is.na(type))) {
        stop("'type' should be one of 'a', 'l', 'p', 'b', 'confint'",
             " or combinations thereof.", call.=FALSE)
    }

    ##- chrTitle
    if (length(chrTitle) != 1) {
        stop("'chrTitle' must be either boolean or character vector of ",
            "length one.", call.=FALSE)
    }

    if (is.na(chrTitle) | is.null(chrTitle)) {
        stop("'chrTitle' must be either boolean or character vector of ",
             "length one.", call.=FALSE)
    }

    if (is.logical(chrTitle)) {
        if (chrTitle) {
            chrTitle <- paste0("Chromosome ",
                                as.character(seqnames(region)))
            showTitles <- TRUE
        } else {
            showTitles <- FALSE
        }
    } else {
        showTitles <- TRUE
    }

    ##- trNames
    if (length(trNames) != 2) {
        stop("'trNames' must be a character vector of length 2.",
             call.=FALSE)
    }

    ##- legend
    if (length(legend) != 1 || !is.logical(legend)) {
        stop("'legend' must be either TRUE or FALSE.", call.=FALSE)
    }

    ##- end checking ---------------------------------------------------------#


    ##- initialization of variables ------------------------------------------#
    ##------------------------------------------------------------------------#

    if (any(c("l", "b") %in% type)) {
        groups <- object@sampleInfo$SampleName
        ord <- order(object@sampleInfo$Condition)
        colGroup <- rep(colGroup, table(object@sampleInfo$Condition))[ord]
    } else {
        groups <- object@sampleInfo$Condition
    }

    chr <- unique(as.vector(seqnames(region)))

    if (length(chr) > 1) {
        warning("multiple chromosomes in region. Only the chromosome in",
                " the first region is taken.", call.=FALSE)
    }

    chr <- chr[1]
    region <- region[seqnames(region) == chr]
    regChr <- rep(chr, length(region))
    id <- c("down", "up")
    id <- id[as.numeric(mcols(region)$log2FoldChange > 0) + 1]

    startReg <- start(region)
    endReg <- end(region)
    strandReg <- strand(region)

    minStart <- min(startReg) - flankReg
    maxEnd <- max(endReg) + flankReg
    len <- maxEnd - minStart + 1

    if (!is.null(annot)) {
        annot <- annot[seqnames(annot) == chr]
        ov <- unique(queryHits(findOverlaps(annot, region)))

        if (length(ov) > 0) {
            annot <- annot[ov]

            if (allSignReg) {
                tmp <- object@regions
                tmp <- tmp[unique(queryHits(findOverlaps(tmp, annot)))]
                idEq <- queryHits(findOverlaps(tmp, region, type = "equal"))
                region <- c(region, tmp[seq_along(tmp)][-idEq])
                regChr <- rep(chr, length(region))
                id <- c("down", "up")
                id <- id[as.numeric(mcols(region)$log2FoldChange > 0) + 1]

                startReg <- start(region)
                endReg <- end(region)
                strandReg <- strand(region)
            }

            minStart <- min(c(start(region), start(annot))) - flankReg
            maxEnd <- max(c(end(region), end(annot))) + flankReg
            len <- maxEnd - minStart + 1

        } else {
            annot <- NULL
        }
    }

    if (length(id) == 0) id <- NULL

    if (normCvg) {
        cvgMat <- as(lapply(cvgNormalization(object), `[[`, chr), "DataFrame")
    } else {
        cvgMat <- as(lapply(coverages(object), `[[`, chr), "DataFrame")
    }

    cvgMat <- cvgMat[minStart:maxEnd, ]

    gr <- GRanges(seqnames = Rle(chr, len),
                  ranges = IRanges(seq(minStart, maxEnd), width = 1),
                  strand = Rle(strand("*"), len),
                  data = cvgMat)


    ##- Tracks ---------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    trackList <- list()

    if (showTitles) {
        pos <- round((startReg + endReg) / 2)
        titleTrack <- AnnotationTrack(start = pos - 1, end = pos + 1,
                                      chromosome = regChr, id = chrTitle,
                                      name = "", fill = "white",
                                      col = "white",
                                      fontcolor.item = "darkgray",
                                      background.title = "white",
                                      featureAnnotation = "id")
        trackList <- c(trackList, titleTrack)
    }

    trackList <- c(trackList, GenomeAxisTrack())

    if (!is.null(annot)) {
        track <- GeneRegionTrack(as.data.frame(annot),
                                 chromosome = chr,
                                 fill = fillAnnot,
                                 name = "")
        trackList <- c(trackList, track)
    }

    track <- AnnotationTrack(start = startReg,
                             end = endReg,
                             strand = strandReg,
                             chromosome = regChr,
                             name = trNames[1],
                             id = id)

    if (!is.null(id)) feature(track) <- id
    trackList <- c(trackList, track)

    track <- DataTrack(gr,
                       groups = groups,
                       legend = legend,
                       col = colGroup,
                       name = trNames[2])
    trackList <- c(trackList, track)


    ##- plot tracks ----------------------------------------------------------#
    ##------------------------------------------------------------------------#
    plotTracks(trackList,
               type = type,
               transcriptAnnotation = featureAnnot,
               shape = "arrow",
               up = fillReg[1],
               down = fillReg[2],
               ...)

    return(invisible(trackList))
}
