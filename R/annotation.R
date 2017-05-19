#' Segmentation using an annotation file.
#'
#' @param fileName The annotation file name in GFF/GTF format.
#' @param source   If not NULL, only lines with this source (2nd field) are imported.
#' @param feature  If not NULL, only lines with this feature (3rd field) are imported.
#' @param name     If not NULL, use this tag as annotation name.
#' @return A GRanges.
#' @examples
#' dir        <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation <- readAnnotation(gtfFile, source="miRNA", feature="gene", name="gene_name")
#'
#' @export
readAnnotation <- function(fileName, source=NULL, feature=NULL, name=NULL) {
    annotation <- import(fileName, feature.type=feature)
    if (! is.null(source)) {
        annotation <- annotation[mcols(annotation)$source == source]
    }
    names(annotation) <- paste("annotation", seq(length(annotation)), sep="_")
    if (! is.null(name)) {
        names(annotation) <- mcols(annotation)[[name]]
    }
    mcols(annotation) <- NULL
    message(paste0(c("... ", length(annotation), " elements found...")))
    message("... Annotation step done.")
    return (annotation)
}

#' Segmentation using an annotation file that contains every genomic feature; it extracts the miRNAs.
#'
#' @param fileName The annotation file name in GFF/GTF format.
#' @return A GRanges.
#' @examples
#' dir        <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gtfFile    <- file.path(dir, "Homo_sapiens.GRCh38.76.gtf")
#' annotation <- readWholeGenomeAnnotation(gtfFile)
#'
#' @export
readWholeGenomeAnnotation <- function(fileName) {
    return(readAnnotation(fileName,
                          source="miRNA", feature="gene", name="gene_name"))
}

#' Segmentation using an miRBase annotation file and use precursor miRNAs.
#'
#' @param fileName The annotation file name in GFF/GTF format.
#' @return A GRanges.
#' @examples
#' dir        <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gffFile    <- file.path(dir, "mirbase21_GRCh38.gff3")
#' annotation <- readMiRBasePreAnnotation(gffFile)
#'
#' @export
readMiRBasePreAnnotation <- function(fileName) {
    return(readAnnotation(fileName,
                          feature="miRNA_primary_transcript", name="Name"))
}

#' Segmentation using an miRBase annotation file and use mature miRNAs.
#'
#' @param fileName The annotation file name in GFF/GTF format.
#' @return A GRanges.
#' @examples
#' dir        <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' gffFile    <- file.path(dir, "mirbase21_GRCh38.gff3")
#' annotation <- readMiRBaseMatureAnnotation(gffFile)
#'
#' @export
readMiRBaseMatureAnnotation <- function(fileName) {
    return(readAnnotation(fileName, feature="miRNA", name="Name"))
}

#' Segmentation using an annotation file.
#'
#' @param object An \code{srnadiff} object.
#' @return A GRanges.
runAllAnnotation <- function(object) {
    if ((is.null(object@annotation)) ||
        (length(object@annotation) == 0) ||
        (object@skipAnnotation)) {
        return(GRanges())
    }
    return (object@annotation)
}
