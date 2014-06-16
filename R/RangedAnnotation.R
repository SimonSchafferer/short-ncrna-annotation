#'@title Ranged Annotation Virtual Class, currently superclass of: GeneAnnotation, RangesAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'  }
#' @name RangedAnnotation-class
#' @export
setClass("RangedAnnotation", representation(annotName = "character", annotationGR = "GRanges", annotationMap="Hits", "VIRTUAL"))

#' Initialisation method for Ranged Annotation class
#'
#' @param .Object GeneAnnotation, GRangesBasedAnnotation
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
			
			return(.Object)
			
		})
#' Constructor method for Ranged Annotation class
#'
#' @param .Object GeneAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
RangedAnnotation = function( annotName, pathToFile, filename, queryGR, ...  ){
  new( "RangedAnnotation", annotName, pathToFile, filename, queryGR, ... )
}


# Accessor Methods --------------------------------------------------------

#' @param RangedAnnotation
#' @return Slot annotationGR
#' @rdname annotName-method
#' @export
setMethod("annotName",signature(object="RangedAnnotation"),function(object) {
  slot(object, "annotName")
})

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

# Other RangedAnnotation Methods ------------------------------------------

#' Import RangedAnnotation
#' @docType methods
#' @rdname importRangedAnnotation-method
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


# GRangesBasedAnnotation class (subclass of RangedAnnotation) --------------------

#'@title GRangesBasedAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'  }
#' @name GRangesBasedAnnotation-class
#' @export
setClass("GRangesBasedAnnotation", contains = "RangedAnnotation")

#' Constructor method for Ranged Annotation class
#'
#' @param .Object GeneAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
GRangesBasedAnnotation = function(  annotName, pathToFile,filename,  queryGR ){
  new( "GRangesBasedAnnotation", annotName, pathToFile, queryGR )
}

#' @param .object GRangesBasedAnnotation
#' @return pathToFile Path to a file supported by rtracklayer/GRanges
#' @rdname importRangedAnnotation-method
setMethod("importRangedAnnotation", "GRangesBasedAnnotation", function(object, pathToFile, filename){
  
  otherAnnot = tryCatch({
    message("...loading library ...")
    loadLibraryLocallyAs(name="otherAnnot", filename=file.path( pathToFile, filename) )
  },warning = function(w){
    message(paste("Probably not an RObject file -> trying to import range based file (bed, gtf, ...)"))
    message(w)
    otherAnnot = rtracklayer::import( file.path( pathToFile, filename), asRangedData = FALSE ) 
    return(otherAnnot)
  },error=function(e){
    message(paste("Probably not an RObject file -> trying to import range based file (bed, gtf, ...)"))
    message(e)
    otherAnnot = rtracklayer::import( file.path( pathToFile, filename), asRangedData = FALSE )
    return(otherAnnot)
  })
  
  if( !is(otherAnnot,"GRanges") ){
    stop("Probably not a GRanges object")
  }
  
  return(otherAnnot)
}  )


# Gene Annotation (subclass of Ranged Annotation) -------------------------

#'@title GeneAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'  }
#' @name GeneAnnotation-class
#' @export
setClass("GeneAnnotation", contains = "RangedAnnotation")

#' @param .object GeneAnnotation
#' @return pathToFile Path to a file supported by rtracklayer/GRanges
#' @rdname importRangedAnnotation-method
setMethod("importRangedAnnotation", "GeneAnnotation", function(object, pathToFile, filename){
  
  geneDF = NULL
  #Try to load r file or read the csv file
  geneDF = tryCatch({
    message("...loading library ...")
    loadLibraryLocallyAs(name="geneDF", filename=file.path(pathToFile,filename) )
  },warning = function(w){
    message("Does not seem to be an RObject")
  },error=function(e){
    message(e)
    return(null)
  })
  
  if(is.null(geneDF) ){
    geneDF = tryCatch({
      message("...creating R object from gtf file ...")
      return( createRObject_gtf( pathToGTF=pathToFile, filename=filename ) )
    },warning = function(w){
      message(w)
    },error=function(e){
      message("Cannot prepare gtf file")
      message(e)
      return(NULL)
    })
  }

  if( is.null(geneDF) ){
    stop("Cannot Read Gene File (should be eather in tab separated format or Robject!)")
  }
  return( geneDF )
}  )



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

#' @param GeneAnnotation
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "GeneAnnotation"), function(object,... ){
    gro = annotationGR(object)
    return(as.data.frame(gro))
})

#' @param GRangesBasedAnnotation
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "GRangesBasedAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})

#' @param GRanges
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "GRanges"), function(object,... ){
  return(as.data.frame(object))
})

