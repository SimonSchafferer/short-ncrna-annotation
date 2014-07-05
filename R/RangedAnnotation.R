#' @title Ranged Annotation Virtual Class, currently superclass of: EnsemblAnnotation, , RangesAnnotation
#' @section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'    \item{\code{slot4}:}{inputGR Input of user that got annotated \code{"GRanges"}}
#'  }
#' @rdname RangedAnnotation-class
#' @export
setClass(Class="RangedAnnotation", representation(annotName = "character", 
                                            annotationGR = "GRanges", 
                                            annotationMap="Hits", 
                                            inputGR="GRanges", "VIRTUAL") )

#' @title Initialisation method for Ranged Annotation class
#'
#' @param .Object EnsemblAnnotation, RefSeqUCSCAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @importFrom rtracklayer import
#' @docType methods
setMethod("initialize", "RangedAnnotation", function(.Object, annotName, pathToFile, filename, queryGR, ... ){
			require(rtracklayer)
			strand(queryGR) = "*"
      message("Importing RangedAnnotation")
			bedToInvestigate <- importRangedAnnotation( .Object, pathToFile, filename, ... )
			
			message("Subsetting to query input")
			bedToInvestigate.sub <- subsetByOverlaps( bedToInvestigate, queryGR )
			#mapping the basicAnnotation file to the subset
			annotationMap <- findOverlaps( queryGR, bedToInvestigate.sub )
			
			if( length(bedToInvestigate.sub) == 0 ){
				bedToInvestigate.sub <- new("GRanges")
				annotationMap <- new("Hits")
			}
			.Object@annotName <-annotName
			.Object@annotationGR <- bedToInvestigate.sub
			.Object@annotationMap <- annotationMap
      .Object@inputGR <- queryGR
			validObject(.Object)
			return(.Object)
			
		})

# Accessor Methods --------------------------------------------------------

#' @title Accessor method annotationGR
#' 
#' @seealso \code{\link{GRanges}}
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @export
#' @docType methods
#' @rdname annotationGR-method
setGeneric("annotationGR", function(object) { standardGeneric("annotationGR") })

#' @rdname annotationGR-method
setMethod("annotationGR",signature(object="RangedAnnotation"),function(object) {
			slot(object, "annotationGR")
		})

#' @title Accessor method annotationMap
#' 
#' @seealso \code{\link{Hits}}
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @export
#' @docType methods
#' @rdname annotationMap-method
setGeneric("annotationMap", function(object) { standardGeneric("annotationMap") })

#' @rdname annotationMap-method
setMethod("annotationMap",signature(object="RangedAnnotation"),function(object) {
			slot(object, "annotationMap")
		})

#' @title Accessor method annotName
#' @export
#' @docType methods
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @rdname annotName-method
setGeneric("annotName", function(object) { standardGeneric("annotName") })

#' @rdname annotName-method
setMethod("annotName",signature(object="RangedAnnotation"),function(object) {
  slot(object, "annotName")
})

#' @title Accessor method inputGR
#' 
#' @seealso \code{\link{GRanges}}
#' @param RangedAnnotation
#' @return Slot inputGR
#' @export
#' @docType methods
#' @rdname inputGR-method
setGeneric("inputGR", function(object) { standardGeneric("inputGR") })

#' @rdname inputGR-method
setMethod("inputGR",signature(object="RangedAnnotation"),function(object) {
  slot(object, "inputGR")
})

#' @title Constructor method for Ranged Annotation class
#'
#' @param .Object EnsemblAnnotation, RefSeqUCSCAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
RangedAnnotation = function( annotName, pathToFile, filename, queryGR, ...  ){
  new( "RangedAnnotation", annotName, pathToFile, filename, queryGR, ... )
}

# ----------------------------------- Other RangedAnnotation Methods ------------------------------------------

#' @title Import Range dAnnotation
#' 
#' @param Object
#' @param pathToFile path to a gtf or rda file
#' @param filename name of the gtf/rda file
#' @rdname importRangedAnnotation-method
#' @docType methods
#' @export
setGeneric("importRangedAnnotation", function( object, pathToFile, filename, ... ) standardGeneric("importRangedAnnotation"))

#' @title Fetch a specific annoation by ID
#' 
#' @docType methods
#' @param RangedAnnotation
#' @return RangedAnnotation object
#' @rdname getAnnotationByID-method
#' @export
setGeneric("getAnnotationByID", function( object, id, ... ) { standardGeneric("getAnnotationByID") })


#' @rdname getAnnotationByID-method
setMethod("getAnnotationByID", signature( "RangedAnnotation"), function(object, id){
			
			annotationMap <- annotationMap(object)
			annotationGR <- annotationGR(object)
			idIndexMap <- which( queryHits(annotationMap) == id )

			if( length(idIndexMap) == 0 ){
				object@annotationGR <- new("GRanges")
				object@annotationMap <- new("Hits")
			} else{
				
				if( length(idIndexMap) == 0){
					object@annotationGR <- new("GRanges")
					object@annotationMap <- new("Hits")
				} else{
					object@annotationGR <- annotationGR[subjectHits(annotationMap)[idIndexMap]]
					object@annotationMap <- annotationMap[idIndexMap]
				}
			}
			return(object)
		})

#' @title Summarize annotation
#' 
#' @param RangedAnnotation class
#' @docType methods
#' @rdname annotationSummary-method
#' @return Annotation data.frames as list
#' @export
setGeneric("annotationSummary", function( object, ... ) { standardGeneric("annotationSummary") })


#' @title Convert Ranges object to data frame
#' 
#' @export
#' @param RangedAnnotation
#' @docType methods
#' @rdname convertRangesToDF-method
#' @return data.frame
setGeneric("convertRangesToDF", function( object, ... ) { standardGeneric("convertRangesToDF") })

#' @rdname convertRangesToDF-method
setMethod("convertRangesToDF", signature("RangedAnnotation"), function(object, ...){
  stop("no convertRangesToDF method for this class")
})

#' @rdname convertRangesToDF-method
setMethod("convertRangesToDF", signature( "GRanges"), function(object,... ){
  return(as.data.frame(object))
})
