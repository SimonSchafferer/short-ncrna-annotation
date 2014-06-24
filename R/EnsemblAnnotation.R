#' @include RangedAnnotation.R

#'@title EnsemblAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'  }
#' @name EnsemblAnnotation-class
#' @export
setClass("EnsemblAnnotation", contains = "RangedAnnotation")


#' Constructor method for Ranged Annotation class
#'
#' @param .Object GeneAnnotation, GRangesBasedAnnotation
#' @param annotName Name of the annotation
#' @param pathToFile file path to either an RObject or specific file type see  \code{"rtracklayer"} of the subject annotation
#' @param queryGR GRanges object of candidates one wants to annotate
#' @export
#' @docType methods
EnsemblAnnotation = function(  annotName, pathToFile, filename, queryGR ){
  new( "EnsemblAnnotation", annotName=annotName, pathToFile=pathToFile, filename, queryGR=queryGR )
}

checkValidityEnsemblAnnotation=function(object) {   
  agr = annotationGR(object)
  
  if( !all(c("type","gene_id","transcript_id", "gene_biotype") %in% colnames(elementMetadata(agr))) ){
    return( "\nEnsembl annotation has to contain a gene_id, transcript_id and gene_biotype column!\n
            gene_id and transcript_id must contain Ensembl Identifier (ENS...)\n
            gene_biotype must contain biotypes e.g. protein_coding, miRNA\n
            type must contain intron, exon, etc." )
  } else{
    return(TRUE)
  }
}
setValidity("EnsemblAnnotation", checkValidityEnsemblAnnotation)

#' importRangedAnnotation
#' @param .object EnsemblAnnotation
#' @param pathToFile path to a gtf file or .rda file
#' @param filename
#' @return Ensembl annotation object rtracklayer/GRanges
#' @rdname importRangedAnnotation-method
setMethod("importRangedAnnotation", "EnsemblAnnotation", function(object, pathToFile, filename){
  
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
      return( createRObject_gtf( pathToGTF=pathToFile, filename=filename, type = "ENSEMBL" ) )
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


#' @param EnsemblAnnotation
#' @return data.frame
#' @rdname convertRangesToDF-method
#' @export
setMethod("convertRangesToDF", signature( "EnsemblAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})
