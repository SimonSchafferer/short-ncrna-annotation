#' @include RangedAnnotation.R

#'@title RefSeqUCSCAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'    \item{\code{slot4}:}{inputGR Input of user that got annotated \code{"GRanges"}}
#'  }
#' @rdname RefSeqUCSCAnnotation-class
#' @docType methods
#' @export
setClass("RefSeqUCSCAnnotation", contains = "RangedAnnotation")


#' @title Constructor method for RefSeqUCSCAnnotation Annotation class
#'
#' @param .Object RefSeqUCSCAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
RefSeqUCSCAnnotation = function(  annotName, pathToFile,filename,  queryGR ){
  new( "RefSeqUCSCAnnotation", annotName=annotName, pathToFile=pathToFile,filename=filename, queryGR=queryGR )
}

checkValidityRefSeqUCSCAnnotation=function(object) {   
  agr = annotationGR(object)
  
  if( !all(c("gene_id","transcript_id") %in% colnames(elementMetadata(agr))) ){
    return( "RefSeqUCSCAnnotation annotation has to contain a gene_id, transcript_id column!\n
            gene_id and transcript_id must contain RefSeq Identifier (NR_ or NM_ etc.)" )
  } else{
    return(TRUE)
  }
}
setValidity("RefSeqUCSCAnnotation", checkValidityRefSeqUCSCAnnotation)

#' @rdname importRangedAnnotation-method
setMethod("importRangedAnnotation",  signature(object="RefSeqUCSCAnnotation"), function(object, pathToFile, filename, ... ){
  
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
      return( createRObject_gtf( pathToGTF=pathToFile, filename=filename, type = "UCSC" ) )
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


#' @rdname convertRangesToDF-method
setMethod("convertRangesToDF", signature( "RefSeqUCSCAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})

#' @rdname annotationSummary-method
setMethod("annotationSummary", signature("RefSeqUCSCAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})
