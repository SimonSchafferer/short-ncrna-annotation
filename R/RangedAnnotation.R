#'@title Ranged Annotation Virtual Class, currently superclass of: EnsemblAnnotation, , RangesAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'  }
#' @name RangedAnnotation-class
#' @export
setClass("RangedAnnotation", representation(annotName = "character", annotationGR = "GRanges", annotationMap="Hits", "VIRTUAL"))
# @include GenericMethodsGlobal.R

#' Initialisation method for Ranged Annotation class
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
			validObject(.Object)
			return(.Object)
			
		})

# Accessor Methods --------------------------------------------------------

#' Accessor method annotationGR
#' 
#' @seealso \code{\link{GRanges}}
#' 
#' @export
#' @docType methods
#' @rdname annotationGR-method
setGeneric("annotationGR", function(object) { standardGeneric("annotationGR") })
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @rdname annotationGR-method
#' @export
setMethod("annotationGR",signature(object="RangedAnnotation"),function(object) {
			slot(object, "annotationGR")
		})

#' Accessor method annotationMap
#' 
#' @seealso \code{\link{Hits}}
#' 
#' @export
#' @docType methods
#' @rdname annotationMap-method
setGeneric("annotationMap", function(object) { standardGeneric("annotationMap") })
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @rdname annotationMap-method
#' @export
setMethod("annotationMap",signature(object="RangedAnnotation"),function(object) {
			slot(object, "annotationMap")
		})

#' Accessor method annotName
#' @export
#' @docType methods
#' @rdname annotName-method
setGeneric("annotName", function(object) { standardGeneric("annotName") })
#' @param RangedAnnotation
#' @return Slot annotationGR
#' @rdname annotName-method
#' @export
setMethod("annotName",signature(object="RangedAnnotation"),function(object) {
  slot(object, "annotName")
})


#' Constructor method for Ranged Annotation class
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

# Other RangedAnnotation Methods ------------------------------------------

#' Import RangedAnnotation
#' @docType methods
#' @rdname importRangedAnnotation-method
#' @export
setGeneric("importRangedAnnotation", function( object, pathToFile, filename, ... ) { standardGeneric("importRangedAnnotation") })

#' Fetch a specific annoation by ID
#' @export
#' @docType methods
#' @rdname getAnnotationByID-method
setGeneric("getAnnotationByID", function( object, id, ... ) { standardGeneric("getAnnotationByID") })
#' @param RangedAnnotation
#' @return RangedAnnotation
#' @rdname getAnnotationByID-method
#' @export
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

# # ' Summarize annotation
# # ' @export
# # ' @docType methods
# # ' @rdname annotationSummary-method
# setGeneric("annotationSummary", function( object, ... ) { standardGeneric("annotationSummary") })

#' Convert Ranges object to data frame
#' @export
#' @docType methods
#' @rdname convertRangesToDF-method
setGeneric("convertRangesToDF", function( object, ... ) { standardGeneric("convertRangesToDF") })

#' @param RangedAnnotation
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature("RangedAnnotation"), function(object, ...){
  stop("no convertRangesToDF method for this class")
})

#' @param GRanges
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "GRanges"), function(object,... ){
  return(as.data.frame(object))
})




