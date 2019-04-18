#' Reads and parses GFF/GTF files
#'
#' \code{readAnnotation} reads and parses content of GFF/GTF files
#' and stores annotated genomic features (regions) in a \code{GRanges}
#' object.
#'
#' \code{feature} and \code{source} can be \code{NULL}. In this case, no
#' selection is performed and all content into the file is imported.
#' If \code{tagName} is \code{NULL}, then a systematic name
#' (\code{annot.N}) is given to elements of the \code{GRanges} object.
#'
#' @param fileName A path, URL or connection to the GFF/GTF annotation
#'                 file. Compressed files (\code{"gz"}, \code{"bz2"} and
#'                 \code{"xz"}) are handled transparently.
#' @param feature  \code{NULL} (the default) or a character vector of valid
#'                 feature types. If not \code{NULL}, then only the features
#'                 of the specified type(s) are imported.
#' @param source   \code{NULL} (the default) or a character vector of valid
#'                 source types. If not \code{NULL}, then only the sources
#'                 of the specified type(s) are imported.
#' @param tagName  \code{NULL} (the default). If not \code{NULL}, use this tag
#'                 as systematic name for the elements of the \code{GRanges}
#'                 object.
#'
#' @return         A \code{GRanges} object with annotated regions information.
#'
#' @seealso
#' \code{\link[rtracklayer]{GFFFile-class}}
#'
#' @examples
#' \dontrun{
#'
#' ##-----------------------------------------------------------------------
#' ## Extraction of miRNAs using an GTF annotation file
#' ##-----------------------------------------------------------------------
#' basedir  <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gtfFile  <- file.path(basedir, "Homo_sapiens.GRCh38.76.gtf.gz")
#' annotReg <- readAnnotation(gtfFile, feature="gene", source="miRNA")
#' annotReg
#'
#' ##-----------------------------------------------------------------------
#' ## Extraction of mature miRNAs using a miRBase-formatted file
#' ##-----------------------------------------------------------------------
#' basedir  <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gffFile  <- file.path(basedir, "mirbase21_GRCh38.gff3")
#' annotReg <- readAnnotation(gffFile, feature="miRNA", tagName="miRNA")
#' annotReg
#'
#' ##-----------------------------------------------------------------------
#' ## Extraction of precursor miRNAs using a miRBase-formatted file
#' ##-----------------------------------------------------------------------
#' basedir  <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gffFile  <- file.path(basedir, "mirbase21_GRCh38.gff3")
#' annotReg <- readAnnotation(gffFile, feature="miRNA_primary_transcript")
#' annotReg
#' }
#'
#' @export
readAnnotation <- function(fileName, feature=NULL, source=NULL, tagName=NULL) {

    annot <- import(fileName)

    ##- checking input arguments and initialize variables --------------------#
    ##------------------------------------------------------------------------#
    ##- feature
    if (!is.null(feature)) {
        idFeature <- NULL
        available <- unique(mcols(annot)$type)
        idTrue <- (feature %in% available)

        if (any(!idTrue)) {
            stop("'feature' should be one or several of:\n",
                paste(available, collapse = ", "), ".", call.=FALSE)
        }

        idFeature <- (mcols(annot)$type %in% feature)
        annot <- annot[idFeature]
    }

    ##- source
    if (!is.null(source)) {
        idSource <- NULL
        available <- unique(mcols(annot)$source)
        idTrue <- (source %in% available)

        if (!any(idTrue)) {
            stop("'source' should be one or several of:\n",
                paste(available, collapse = ", "), ".", call.=FALSE)
        }

        idSource <- (mcols(annot)$source %in% source)
        annot <- annot[idSource]
    }

    ##- tagName
    if (!is.null(tagName)) {
        if (!isSingleString(tagName)) {
            stop("'tagName' should be a single character.", call.=FALSE)
        }
        names(annot) <- paste(tagName, 1:length(annot), sep=".")
    } else {
        names(annot) <- paste("annot", 1:length(annot), sep=".")
    }

    ##- end checking ---------------------------------------------------------#

    idAllNAs <- apply(is.na(as.data.frame(mcols(annot))), 2, all)
    mcols(annot) <- mcols(annot)[, !idAllNAs]

    message(paste0(c("... ", length(annot), " elements found...")))
    message("... Annotation step done.")

    return(annot)
}
