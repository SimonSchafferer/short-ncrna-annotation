#'@title SequenceBasedAnnotation Virtual Class, currently superclass of: ncbiBlastAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationL Overlap with annotation \code{"list"}}
#'  }
#' @name SequenceBasedAnnotation-class
#' @export
setClass("SequenceBasedAnnotation", representation(annotName = "character", annotationL = "list", "VIRTUAL"))
#' Initialisation method for SequenceBasedAnnotation class
#'
#' @param .Object ncbiBlastAnnotation
#' @param annotName of the annotation file
#' @param annotationL list of annotations
#' @export
#' @docType methods
#' 
setMethod("initialize", "SequenceBasedAnnotation", function(.Object, annotName, pathToFile, inputStringSet, ... ){
  require(Biostrings)
  annotationL = executeAnnotation(.Object, pathToFile, inputStringSet, ... )
  .Object@annotName = annotName
  .Object@annotationL = annotationL
  return(.Object)
})


#' @param SequenceBasedAnnotation
#' @return Slot annotName
#' @rdname annotName-method
#' @export
setMethod("annotName",signature(object="SequenceBasedAnnotation"),function(object) {
  slot(object, "annotName")
})

#' Accessor method annotationL
#' @export
#' @docType methods
#' @rdname annotationL-method
setGeneric("annotationL", function(object) { standardGeneric("annotationL") })
#' @param SequenceBasedAnnotation
#' @return Slot annotationL
#' @rdname annotationL-method
#' @export
setMethod("annotationL",signature(object="SequenceBasedAnnotation"),function(object) {
  slot(object, "annotationL")
})

#' executeAnnotation
#' 
#' Executes the annotation based on a XStringSet object
#' 
#' @param object currently supported: NcbiBlastAnnotation
#' @param pathToFile path to an annotation File (currently ncbi makeblast db formatted database)
#' @param inputStringSet XStringSet of the sequences one wants to annotate
#' @param ... special options for subclass defintions of executeAnnotation
#' @docType methods
#' @rdname executeAnnotation-method
setGeneric("executeAnnotation", function( object, pathToFile, inputStringSet, ... ) { standardGeneric("executeAnnotation") })

# ncbiBlastAnnotation -----------------------------------------------------

#'@title NcbiBlastAnnotation Class
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationL Overlap with annotation \code{"list"}}
#'  }
#' @name NcbiBlastAnnotation-class
#' @export
setClass("NcbiBlastAnnotation", contains = "SequenceBasedAnnotation")
#' executeAnnotation
#' 
#' Executes the annotation based on a XStringSet object
#' 
#' @param object currently supported: NcbiBlastAnnotation
#' @param pathToFile path to an annotation File (currently ncbi makeblast db formatted database)
#' @param inputStringSet XStringSet of the sequences one wants to annotate
#' @param ncbiBlast specific parameters see: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#' @docType methods
#' @rdname executeAnnotation-method
#' @import XML
#' @import Biostrings
setMethod("executeAnnotation", "NcbiBlastAnnotation", function(object,pathToFile, inputStringSet, word_size=11, num_threads=2, dust=TRUE ){
  require(XML)
  require(Biostrings)
  max_target_seqs = 1 #Currently not modifiable, taking best hit!
  bestHit = 1 #Currently not modifiable, taking best hit!
  #Nucleotide-Nucleotide BLAST 2.2.25+
  tmpFile = tempfile()
  message( paste0("Writing input sequences temporarily to ",tmpFile) )
  Biostrings::writeXStringSet(x=inputStringSet, filepath=tmpFile )
  dust = if(dust){"yes"}else{"no"}
  
  blastcmd <- paste("blastn -query", tmpFile ,"-db", pathToFile ,"-max_target_seqs",
                    max_target_seqs,"-word_size",word_size,"-num_threads",num_threads,"-outfmt \"5\"", " -dust ", dust)
  
  message(paste0( "Executing Blast command: ", blastcmd ))
  
  blastout <- system(blastcmd, intern = TRUE)
  blastoutxml = XML::xmlParse(blastout)
  blastoutxml.iter <- XML::xpathSApply(blastoutxml, "//Iteration"); 
  blastdfname <- c( )
  xmlDocListHit <- lapply( blastoutxml.iter, function(x){
    message("... parsing the XML output ...")
    x <- xmlDoc(x)
    blastdfname <<- c(blastdfname,xmlToDataFrame(xpathSApply( x, "//Iteration" ))[,c("Iteration_query-ID")])  					
    blastdf <- xmlToDataFrame(xpathSApply(x, "//Hit"))
    
    if(dim( blastdf )[1] != 0){
      #Hit found -> else the empty data frame is returned
      blastdf <- blastdf[,c("Hit_num","Hit_id","Hit_def", "Hit_accession", "Hit_len")]
      tmpblastdf <- xmlToDataFrame(xpathSApply(x, "//Hsp"))[,c("Hsp_num","Hsp_bit-score", "Hsp_score", "Hsp_evalue", "Hsp_query-from", "Hsp_query-to", "Hsp_hit-from", "Hsp_hit-to", "Hsp_gaps", "Hsp_align-len", "Hsp_qseq", "Hsp_hseq")]
      #Currently delete the second annotation in one sequence -> automatically ranked by bit-score
      tmpblastdf <- tmpblastdf[tmpblastdf$Hsp_num == bestHit,]
      blastdf <- cbind( blastdf, tmpblastdf)
    }
    return( blastdf )
  } )
  names(xmlDocListHit) <- blastdfname
  tryCatch({
    message( paste0("Deleting sequences from ",tmpFile) )
    file.remove(tmpFile)
  }, error = function(e){
    message(paste0("File could not be deleted!", e))
  })
  return(xmlDocListHit)
  
})



