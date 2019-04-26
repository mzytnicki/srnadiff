##- bamFiles
setGeneric(name="bamFiles",
            def=function(object) {
                standardGeneric("bamFiles")
            }
)

##- sampleInfo
setGeneric(name="sampleInfo",
            def=function(object) {
                standardGeneric("sampleInfo")
            }
)

##- annotReg
setGeneric(name="annotReg",
            def=function(object) {
                standardGeneric("annotReg")
            }
)

setGeneric(name="annotReg<-",
           def=function(object, value) {
               standardGeneric("annotReg<-")
           }
)

##- chromosomeSizes
setGeneric(name="chromosomeSizes",
            def=function(object) {
                standardGeneric("chromosomeSizes")
            }
)

##- normFactors
setGeneric(name="normFactors",
            def=function(object) {
                standardGeneric("normFactors")
            }
)

setGeneric(name="normFactors<-",
            def=function(object, value) {
                standardGeneric("normFactors<-")
            }
)

##- coverages
setGeneric(name="coverages",
            def=function(object) {
                standardGeneric("coverages")
            }
)

##- regions
setGeneric(name="regions",
            def=function(object, pvalue=1) {
                standardGeneric("regions")
            }
)

##- parameters
setGeneric(name="parameters",
            def=function(object) {
                standardGeneric("parameters")
            }
)

setGeneric(name="parameters<-",
            def=function(object, value) {
                standardGeneric("parameters<-")
            }
)

##- countMatrix
setGeneric(name="countMatrix",
            def=function(object) {
                standardGeneric("countMatrix")
            }
)


##- print.srnadiff_par
setGeneric(name="print.srnadiff_par",
            def=function(object) {
                standardGeneric("print.srnadiff_par")
            }
)
