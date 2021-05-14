##- bamFiles -----------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'bamFiles' slot of an srnadiffExp object
#'
#' The \code{bamFiles} slot holds the full paths to the BAM files as a
#' \code{\link[Rsamtools]{BamFileList}}.
#'
#' @docType methods
#' @name bamFiles
#' @rdname bamFiles
#' @aliases bamFiles bamFiles,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @return The paths to the BAM files.
#' @examples
#' require(Rsamtools)
#'
#' srnaExp <- srnadiffExample()
#' path(bamFiles(srnaExp))
#'
#' @export
setMethod(f="bamFiles", signature="srnadiffExp",
          definition=function(object) {
              object@bamFiles
          }
)


##- sampleInfo ---------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'sampleInfo' slot of an srnadiffExp object
#'
#' The \code{sampleInfo} slot holds the sample information as a
#' \code{data.frame} with three columns labelled \code{FileName},
#' \code{SampleName} and \code{Condition}. The first column is the
#' BAM file name (without extension), the second column the sample
#' name, and the third column the condition to which sample belongs.
#' Each row describes one sample.
#'
#' @docType methods
#' @name sampleInfo
#' @rdname sampleInfo
#' @aliases sampleInfo sampleInfo,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @return A table containing information on the samples
#' @examples
#' srnaExp <- srnadiffExample()
#' sampleInfo(srnaExp)
#'
#' @export
setMethod(f="sampleInfo", signature="srnadiffExp",
          definition=function(object) {
              object@sampleInfo
          }
)


##- annotReg -----------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'annotReg' slot of an srnadiffExp object
#'
#' The \code{annotReg} slot holds the annotated regions as a \code{GRanges}
#' object.
#'
#' @docType methods
#' @name annotReg
#' @rdname annotReg
#' @aliases annotReg annotReg,srnadiffExp-method
#' annotReg<- annotReg<-,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @param value Annotated regions as a \code{GRanges} object.
#' @return The regions given by the user, as a \code{GRanges}.
#' @examples
#' srnaExp <- srnadiffExample()
#' annotReg(srnaExp)
#'
#' @export
setMethod(f="annotReg", signature="srnadiffExp",
            definition=function(object) {
                if (is.null(object@annotReg)) {
                    cat("No 'annotReg' slot found in the srnadiffExp object.")
                } else {
                    object@annotReg
                }
            }
)


#' @name annotReg
#' @rdname annotReg
#' @exportMethod "annotReg<-"
setReplaceMethod("annotReg",
                signature(object="srnadiffExp", value="ANY"),
                    function(object, value) {

                        ##- checking input value -----------------------------#
                        ##----------------------------------------------------#
                        if (!is(value, "GRanges")) {
                            stop("'value' must be an object of class",
                                "'GRanges'.", call.=FALSE)
                        }
                        ##- end check ----------------------------------------#

                        object@annotReg <- value
                        object
                    }
)


##- chromosomeSizes ----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'chromosomeSizes' slot of an srnadiffExp object
#'
#' The \code{chromosomeSizes} slot holds the sizes of the chromosomes
#' as a named vector with chromosome names.
#'
#' @docType methods
#' @name chromosomeSizes
#' @rdname chromosomeSizes
#' @aliases chromosomeSizes chromosomeSizes,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @return A list containing the chromosome sizes.
#' @examples
#' srnaExp <- srnadiffExample()
#' chromosomeSizes(srnaExp)
#'
#' @export
setMethod(f="chromosomeSizes", signature="srnadiffExp",
          definition=function(object) {
              object@chromosomeSizes
          }
)


##- normFactors --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'normFactors' slot of an srnadiffExp object
#'
#' The \code{normFactors} slot holds the normalization factors as a
#' named vector with sample names.
#'
#' The \code{normFactors} vector assigns to each sample coverage a value,
#' the normalization factor, such that count values in each sample coverage
#' can be brought to a common scale by dividing by the corresponding
#' normalization factor. This step is also called normalization, its purpose
#' is to render coverages (counts) from different samples, which may have
#' been sequenced to different depths, comparable. Normalization factors
#' are estimated using the \emph{median ratio method} described by
#' Equation 5 in Anders and Huber (2010). Alternative normalization factor
#' estimators can also be supplied using the assignment function
#' \code{sizeFactors<-}.
#'
#' @references
#' Simon Anders, and Wolfgang Huber (2010). Differential expression analysis
#' for sequence count data. \emph{Genome Biology}, 11:106.
#' @docType methods
#' @name normFactors
#' @rdname normFactors
#' @aliases normFactors normFactors,srnadiffExp-method
#' normFactors<- normFactors<-,srnadiffExp,numeric-method
#' @param object An \code{srnadiffExp} object.
#' @param value A numeric vector, one size factor for each sample
#'              in the coverage data.
#' @return The normalization factors, in a list.
#' @examples
#' srnaExp <- srnadiffExample()
#' normFactors(srnaExp)
#'
#' @export
setMethod(f="normFactors", signature="srnadiffExp",
            definition=function(object) {
                if (is.null(object@normFactors)) {
                    message("No 'normFactors' slot found in the srnadiffExp",
                            " object. Run srnadiff first or assign a vector",
                            " of normalization factors.")
                    } else {
                        object@normFactors
                    }
            }
)

#' @name normFactors
#' @rdname normFactors
#' @exportMethod "normFactors<-"
setReplaceMethod("normFactors",
                signature(object="srnadiffExp", value="numeric"),
                    function(object, value) {

                        ##- checking input value -----------------------------#
                        ##----------------------------------------------------#
                        n <- nrow(sampleInfo(object))

                        if (length(value) != n) {
                            stop("'value' must be a vector of length ", n, ".",
                                call.=FALSE)
                        }

                        if (any(is.na(value))) {
                            stop("NA values in 'value'.", call.=FALSE)
                        }

                        if (!all(is.finite(value))) {
                            stop("Infinite values in 'value'.", call.=FALSE)
                        }

                        if (!all(value > 0)) {
                            stop("'value' must be a vector of positive",
                                " values.", call.=FALSE)
                        }

                        factorNames <- names(value)
                        sampleName <- sampleInfo(object)$SampleName

                        if (is.null(factorNames)) {
                            names(value) <- sampleInfo(object)$SampleName
                        } else {
                            if (any(factorNames != sampleName)) {
                                stop("'value' names must be equal to",
                                    " 'SampleName' column in sampleInfo.",
                                    call.=FALSE)
                            }
                        }
                        ##- end check ----------------------------------------#

                        object@normFactors <- value
                        object
                    }
)


##- coverages ----------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'coverages' slot of an srnadiffExp object
#'
#' The \code{coverages} slot holds the coverages as a named
#' \code{\link[IRanges]{RleList}} object with one coverage vector
#' per seqlevel in the set of genomic alignments from the BAM files.
#'
#' @docType methods
#' @name coverages
#' @rdname coverages
#' @aliases coverages coverages,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @return The coverages, as a list of \code{RleList}.
#' @examples
#' srnaExp <- srnadiffExample()
#' coverages(srnaExp)
#'
#' @export
setMethod(f="coverages", signature="srnadiffExp",
          definition=function(object) {
              object@coverages
          }
)


##- regions ------------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Extracts differentially expressed regions of an srnadiffExp object
#'
#' This function extracts the differentially expressed regions from
#' \code{\link{srnadiffExp}} object ranked by p-value.
#'
#' @docType methods
#' @name regions
#' @rdname regions
#'
#' @aliases regions regions,srnadiffExp-method
#'
#' @param object An \code{srnadiffExp} object.
#' @param pvalue Numeric cutoff value for adjusted p-values. Only regions with
#'               adjusted p-values equal or lower than specified are returned.
#'               Default to 1, all regions are returned.
#'
#' @return A \code{GenomicRanges} object of the selected differentially
#' expressed regions.
#'
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#' regions(srnaExp, pvalue = 0.05)
#'
#' @export
setMethod(f="regions", signature="srnadiffExp",
            definition=function(object, pvalue=1) {
                if (is.null(object@regions)) {
                    message("No 'regions' slot found in the srnadiffExp",
                            " object. Run srnadiff first.")
                } else {
                    ##- checking input value ---------------------------------#
                    ##--------------------------------------------------------#
                    if (length(pvalue) != 1) {
                        stop("'pvalue' must be a single value.", call.=FALSE)
                    }

                    if (is.null(pvalue) || !is.numeric(pvalue) ||
                        !is.finite(pvalue)) {
                        stop("'pvalue' value must be numeric.", call.=FALSE)
                    }

                    if ((pvalue > 1) || (pvalue < 0)) {
                        stop("'pvalue' value ", pvalue, ", outside the",
                            " interval [0,1].", call.=FALSE)
                    }
                    ##- end check -------------------------------------------#

                    reg <- object@regions[mcols(object@regions)$padj <= pvalue]
                    ord <- order(mcols(reg)$padj)
                    reg <- reg[ord]
                    return(reg)
                }
            }
)


##- parameters ---------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'parameters' slot of an srnadiffExp object
#'
#' The \code{parameters} slot holds the parameter values
#' used in an sRNA-diff approach as a named \code{list}. Default values
#' exist for parameters, but these can also be supplied as input
#' values in the \code{useParameters} argument of the \code{\link{srnadiff}}
#' function or using the assignment function \code{\link{parameters<-}}.
#'
#' Parameters in an sRNA-diff approach.
#'
#' \subsection{Global parameters}{
#'    \describe{
#'        \item{\code{minDepth}}{The cutoff to filter the base-level
#'              coverage. Bases where at least one sample has (normalized)
#'              coverage greater than \code{minDepth} be been retained.
#'              Default to \code{10}.}
#'        \item{\code{minSize}}{The minimum size (in base-pairs) of the
#'              regions to be found. Default to \code{18}.}
#'        \item{\code{maxSize}}{The maximum size (in base-pairs) of the
#'              regions to be found. Default to \code{1000000}.}
#'        \item{\code{minGap}}{The minimum gap between regions. Regions
#'              separated by a gap of at most \code{minGap} positions
#'              are merged. Default to \code{100}.}
#'        \item{\code{maxDiff}}{The maximum number of different bases between
#'              two regions. Near-identical regions are collapsed.
#'              Only regions with at most \code{maxDiff} different
#'              positions are considered identicals and are collapsed
#'              into one single region. Default to \code{20}.}
#'        \item{\code{minOverlap}}{This parameters is used in the construction
#'              of the \code{\link{countMatrix}} matrix. Only reads (ranges)
#'              with a minimum of \code{minOverlap} overlapping each expressed
#'              region are considered to be overlapping. Default to \code{10}.}
#'    }
#' }
#'
#' \subsection{Parameters for the HMM method}{
#'    \describe{
#'        \item{\code{noDiffToDiff}}{Initial transition probability from
#'              no differentially expressed state to differentially expressed.
#'              Default to \code{0.001}.}
#'        \item{\code{diffToNoDiff}}{Initial transition probability from
#'              differentially expressed state to no differentially expressed.
#'              Default to \code{0.000001}.}
#'        \item{\code{emission}}{Emission probability. Default to \code{0.9}.}
#'        \item{\code{emissionThreshold}}{Emission threshold. A real number
#'              between \code{0} and \code{1}. Default to \code{0.1}.}
#'    }
#' }
#'
#' \subsection{Parameters for the Naive and IR methods}{
#'    \describe{
#'        \item{\code{minLogFC}}{The minimun sliding threshold used in the
#'              naive and IR method. Default to \code{0.5}.}
#'    }
#' }
#'
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,srnadiffExp-method
#' parameters<- parameters<-,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @param value  A named \code{list} containing valid parameters. See details.
#' @return The named list of the parameters used in the analysis.
#' @seealso
#' \code{useParameters} argument in \code{\link{srnadiff}} function.
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#' print(parameters(srnaExp))
#'
#' parameters(srnaExp) <- list("minSize" = 1, "maxSize" = 1500)
#'
#' @export
setMethod(f="parameters", signature="srnadiffExp",
            definition=function(object) {
                if (is.null(object@parameters)) {
                    message("No 'parameters' slot found in the srnadiffExp",
                            " object. Run srnadiff first or assign a named",
                            " list of valid parameters. See help(parameters)",
                            " for details.")
                } else {
                    object@parameters
                    class(object@parameters) <- "srnadiff_par"
                    return(invisible(object@parameters))
                }
            }
)

#' @name parameters
#' @rdname parameters
#' @exportMethod "parameters<-"
setReplaceMethod("parameters",
                signature(object="srnadiffExp", value="ANY"),
                function(object, value) {

                    ##- checking input value ---------------------------------#
                    ##--------------------------------------------------------#
                     defaultParNames <- names(srnadiffDefaultParameters)

                     if (!is.null(object@parameters)) {
                         srnadiffDefaultParameters <- object@parameters
                     }

                     if (!is(value, "list")) {
                         print(value)
                         print(typeof(value))
                         print(class(value))
                         stop("'value' must be a named list. See",
                             " help(parameters) for details.", call.=FALSE)
                     }

                     valueNames <- names(value)

                     if (any(duplicated(valueNames))) {
                         stop("duplicate name parameters in 'value'. See",
                              " help(parameters) for details.", call.=FALSE)
                     }

                     if (!all(valueNames %in% defaultParNames)) {
                         stop("'value' must be a named list of valid",
                            " parameters. See help(parameters) for details.",
                            call.=FALSE)
                     }

                     ##- individual parameters
                     srnadiffDefaultParameters[valueNames] <- value
                     checkParameters(srnadiffDefaultParameters)

                     ##- end check -------------------------------------------#

                     object@parameters <- srnadiffDefaultParameters
                     object
                 }
)


##- countMatrix --------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Accessors for the 'countMatrix' slot of an srnadiffExp object
#'
#' The \code{countMatrix} slot holds the count matrix from the DERs found.
#'
#' In the last step of an sRNA-diff approach, the potential DERs is assessed.
#' Reads (including fractions of reads) that overlap each expressed region
#' are counted to arrive at a count matrix with one row per region and one
#' column per sample. This matrix to been used for quantify the statistic
#' signification of the finded regions.
#'
#' @docType methods
#' @name countMatrix
#' @rdname countMatrix
#' @aliases countMatrix countMatrix,srnadiffExp-method
#' @param object An \code{srnadiffExp} object.
#' @return A matrix with the number of reads for each regions, and each sample.
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#' countMatrix(srnaExp)
#'
#' @export
setMethod(f="countMatrix", signature="srnadiffExp",
            definition=function(object) {
                if (is.null(object@countMatrix)) {
                    message("No 'countMatrix' slot found in the srnadiffExp",
                            " object. Run srnadiff first.")
                } else {
                    object@countMatrix
                }
            }
)


##- show ---------------------------------------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname srnadiffExp
#' @param object An \code{srnadiffExp} object.
#' @return The \code{show} method informatively display object contents.
#' @export
setMethod(f="show", signature ="srnadiffExp",
            definition=function(object) {
                cat("Object of class srnadiffExp.\n",
                    "Sample information\n")
                print(object@sampleInfo)
            }
)


##- print method for parameters ----------------------------------------------#
##----------------------------------------------------------------------------#
#' Dispatch print method for the parameters used by an \code{srnadiff} object.
#'
#' @docType methods
#' @name parameters
#' @rdname parameters
#' @aliases parameters parameters,srnadiffExp-method
#' @param x The first element of the parameters used by an \code{srnadiff}
#'          object
#' @param ... The other elements of the parameters
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#' print(parameters(srnaExp))
#'
#' @export
print.srnadiff_par <- function(x, ...) {

    cat("\n Global parameters: \n",
        "------------------ \n")
    df <- data.frame(value = unlist(x[seq_len(6)]))
    print(df)

    cat("\n HMM method parameters: \n",
        "---------------------- \n")
    df <- data.frame(value = unlist(x[7:10]))
    print(df)

    cat("\n IR method parameter: \n",
        "----------------------- \n")
    df <- data.frame(value = unlist(x[12]))
    print(df)

    cat("\n Naive method parameter: \n",
        "----------------------- \n")
    df <- data.frame(value = unlist(x[11]))
    print(df)
}
