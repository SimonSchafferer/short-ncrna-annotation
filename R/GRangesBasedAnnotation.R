#' @include RangedAnnotation.R

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
#' @param .Object EnsemblAnnotation, RefSeqGUCSCAnnoation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
GRangesBasedAnnotation = function(  annotName, pathToFile,filename,  queryGR ){
  new( "GRangesBasedAnnotation", annotName=annotName, pathToFile=pathToFile,filename=filename, queryGR=queryGR )
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

#' @param GRangesBasedAnnotation
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "GRangesBasedAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})
