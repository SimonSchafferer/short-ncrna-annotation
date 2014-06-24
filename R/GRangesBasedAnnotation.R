#' @include RangedAnnotation.R

# GRangesBasedAnnotation class (subclass of RangedAnnotation) --------------------
#'@title GRangesBasedAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'    \item{\code{slot4}:}{inputGR Input of user that got annotated \code{"GRanges"}}
#'  }
#' @rdname GRangesBasedAnnotation-class
#' @docType class
#' @export
setClass("GRangesBasedAnnotation", contains = "RangedAnnotation")

#' @title Constructor method for Ranged Annotation class
#'
#' @param .Object EnsemblAnnotation, RefSeqGUCSCAnnoation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
GRangesBasedAnnotation = function(  annotName, pathToFile,filename,  queryGR ){
  new( "GRangesBasedAnnotation", annotName=annotName, pathToFile=pathToFile,filename=filename, queryGR=queryGR )
}

#' @rdname importRangedAnnotation-method
setMethod("importRangedAnnotation", signature(object="GRangesBasedAnnotation"), function(object, pathToFile, filename, ...){
  
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


#' @rdname convertRangesToDF-method
setMethod("convertRangesToDF", signature( "GRangesBasedAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})


#' @rdname annotationSummary-method
setMethod("annotationSummary", signature("GRangesBasedAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})

