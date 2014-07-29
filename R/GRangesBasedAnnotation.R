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



#' @title feature_Annotation_gr helper function
#' 
#' @param annotGR
#' @param annotMap
#' @param inputGR
#' @param inputIDs
#' @docType methods
#' @rdname feature_annotation_Function-method
#' @return data.frame
feature_annotation_Function_gr = function(annotGR, annotMap, inputGR, inputIDs){
  
  if( length(annotGR) == 0 ){
    return(data.frame())
  }
  
  alignmentCoverage_featureDF = calculateAlignmentCoverageTwoGRanges( qh=inputGR, sh=annotGR )
  colnames(alignmentCoverage_featureDF) = c("feature_queryCoverage","feature_subjectCoverage")
  
  featureDF = data.frame("InputID"=inputGR$tmpID, stringsAsFactors=FALSE)
  featureDF = cbind( featureDF, as.data.frame( elementMetadata(annotGR) ) )
  featureDF = cbind(featureDF, alignmentCoverage_featureDF)
  featureDF$featureStrand = as.character(strand(annotGR))
  
  missingIDs = inputIDs[which(!inputIDs %in% featureDF$InputID)]
  
  if( length(missingIDs) != 0 ){
    missingDF = data.frame(matrix(nrow = length(missingIDs), ncol = dim(featureDF)[2]), stringsAsFactors=FALSE)
    colnames(missingDF) = colnames(featureDF)
    missingDF$InputID = missingIDs
    featureDF = rbind( featureDF, missingDF)
  }
  rownames(featureDF) = 1:dim(featureDF)[1]
  featureDF = featureDF[order(featureDF$InputID),]
  return(featureDF)
}

#' @rdname getFlatTable-method
setMethod("getFlatTable", signature("GRangesBasedAnnotation"), function(object, ...){
  
  #First defining the genome location and the rest are feature overlaps
  annotgr = annotationGR(object)
  mapgro = annotationMap(object)
  ingr = inputGR(object)
  ingr$tmpID = 1:length(ingr)
  validObject(object)
  
  tmpIngr = ingr[queryHits(mapgro)]
  tmpAnnotGR = annotgr[subjectHits(mapgro)]
  
  #Indices are refering to the original file mapgro
  feature_annotation_map = mapgro
  feature_AnnotGR = tmpAnnotGR
  featureInGR = ingr[queryHits(feature_annotation_map)]
  
  ###############################################
  #  Calling the feature annotation Function
  ###############################################
  featureDF = feature_annotation_Function_gr( annotGR=feature_AnnotGR, 
                                              annotMap=feature_annotation_map, 
                                              inputGR=featureInGR, 
                                              inputIDs=ingr$tmpID)
  return(featureDF)
})

#' @rdname annotationSummary-method
setMethod("annotationSummary", signature("GRangesBasedAnnotation"), function(object,... ){

  featureDF = getFlatTable(object)
  
  rowCoverage = rowSums( featureDF[,c("feature_queryCoverage", "feature_subjectCoverage")] )# 2 is the highest number in rowSums
  rowCoverage = max(rowCoverage,na.rm=TRUE)-rowCoverage
  featureDF$rowCoverage = rowCoverage
  featureDF = featureDF[order(featureDF$InputID, featureDF$rowCoverage),]
  
  featid = featureDF$InputID
  featgt = featureDF$InputID
  featureDFL = split(featgt, featid)
  featureDF_summary = featureDF[ which( !duplicated(featureDF$InputID)),]#always first entry
  featureDF_summary$NumberOfTranscripts =  unlist( lapply( featureDFL, function(x){ 
    return(length(x))
  }))
  
  rownames(featureDF_summary) = 1:dim(featureDF_summary)[1]
  
  return(featureDF_summary)
})

