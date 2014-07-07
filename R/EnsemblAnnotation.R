#' @include GRangesBasedAnnotation.R

#'@title EnsemblAnnotation extends RangedAnnotation
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{annotName of the annotation \code{"character"}}
#'    \item{\code{slot2}:}{annotationGR Overlap with annotation \code{"GRanges"}}
#'    \item{\code{slot3}:}{annotationMap Map between annotation entry and input Entry \code{"Hits"}}
#'    \item{\code{slot4}:}{inputGR Input of user that got annotated \code{"GRanges"}}
#'  }
#' @rdname EnsemblAnnotation-class
#' @aliases EnsemblAnnotation EnsemblAnnotation-class
#' @export
setClass("EnsemblAnnotation", contains = "RangedAnnotation")

#' @title Constructor method for Ranged Annotation class
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

#' @rdname importRangedAnnotation-method
#' @aliases importRangedAnnotation, EnsemblAnnotation
setMethod("importRangedAnnotation", signature(object="EnsemblAnnotation"), function(object, pathToFile, filename, ...){
  
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


#' @rdname convertRangesToDF-method
setMethod("convertRangesToDF", signature( "EnsemblAnnotation"), function(object,... ){
  gro = annotationGR(object)
  return(as.data.frame(gro))
})

#' @title protein_coding_Annotation_ensembl helper function
#' 
#' @param annotGR
#' @param annotMap
#' @param inputGR
#' @param inputIDs
#' @docType methods
#' @rdname protein_coding_AnnotationFunction-method
#' @return data.frame
protein_coding_AnnotationFunction_ensembl = function( annotGR, annotMap, inputGR, inputIDs){
  
  #Subsetting to necessary information
  geneDF = data.frame( 
    "InputID"=inputGR$tmpID, 
    "gene_biotype"=annotGR$gene_biotype,
    "gene_strand"=as.character(strand(annotGR)),
    "gene_type"=annotGR$type,
    "gene_id"=annotGR$gene_id,
    "transcript_id"=annotGR$transcript_id,
    "exonIntron_number"=annotGR$exon_number,
    "gene_name"=annotGR$gene_name,
    "gene_transcript_name"=annotGR$transcript_name, 
    stringsAsFactors=FALSE
  )
  
  missingIDs = inputIDs[which(!inputIDs %in% geneDF$InputID)]
  
  if( length(missingIDs) != 0 ){
    missingDF = data.frame( 
      "InputID"=missingIDs, 
      "gene_biotype"=NA,
      "gene_strand"=NA,
      "gene_type"="intergenic",
      "gene_id"=NA,
      "transcript_id"=NA,
      "exonIntron_number"=NA,
      "gene_name"=NA,
      "gene_transcript_name"=NA, 
      stringsAsFactors=FALSE
    )
    
    geneDF = rbind( geneDF, missingDF)
  }
  rownames(geneDF) = 1:dim(geneDF)[1]
  geneDF = geneDF[order(geneDF$InputID),]
  
  return(geneDF)
}

#' @title feature_Annotation_ensembl helper function
#' 
#' @param annotGR
#' @param annotMap
#' @param inputGR
#' @param inputIDs
#' @docType methods
#' @rdname feature_annotation_Function-method
#' @return data.frame
#' @export
feature_annotation_Function_ensembl = function(annotGR, annotMap, inputGR, inputIDs){
  
  alignmentCoverage_featureDF = calculateAlignmentCoverageTwoGRanges( qh=inputGR, sh=annotGR )
  colnames(alignmentCoverage_featureDF) = c("feature_queryCoverage","feature_subjectCoverage")
  
  #Subsetting to necessary information
  featureDF = data.frame( 
    "InputID"=inputGR$tmpID, 
    "feature_biotype"=annotGR$gene_biotype,
    "feature_strand" = as.character(strand(annotGR)),
    "feature_type"=annotGR$type,
    "feature_id"=annotGR$gene_id,
    "feature_transcript_id"=annotGR$transcript_id,
    "feature_exonIntron_number"=annotGR$exon_number,
    "feature_name"=annotGR$gene_name, 
    stringsAsFactors=FALSE)
  featureDF = cbind(featureDF, alignmentCoverage_featureDF)
  
  missingIDs = inputIDs[which(!inputIDs %in% featureDF$InputID)]
  
  if( length(missingIDs) != 0 ){
    missingDF = data.frame( 
      "InputID"=missingIDs, 
      "feature_biotype"=NA,
      "feature_strand" = NA,
      "feature_type"=NA,
      "feature_id"=NA,
      "feature_transcript_id"=NA,
      "feature_exonIntron_number"=NA,
      "feature_name"=NA,
      "feature_queryCoverage"=NA,
      "feature_subjectCoverage"=NA,
      stringsAsFactors=FALSE
    )
    featureDF = rbind( featureDF, missingDF)
  }
  rownames(featureDF) = 1:dim(featureDF)[1]
  featureDF = featureDF[order(featureDF$InputID),]
  
  return(featureDF)
  
}

#' @rdname getFlatTable-method
setMethod("getFlatTable", signature("EnsemblAnnotation"), function(object, ...){
  
  #First defining the genome location and the rest are feature overlaps
  annotgr = annotationGR(object)
  mapgro = annotationMap(object)
  ingr = inputGR(object)
  ingr$tmpID = 1:length(ingr)
  validObject(object)
  
  tmpIngr = ingr[queryHits(mapgro)]
  tmpAnnotGR = annotgr[subjectHits(mapgro)]
  
  ###############################################
  #  separating the protein-coding gene annotation from the feature annotation
  ###############################################
  
  #Indices are refering to the original file mapgro
  protein_codingIdx = which(tmpAnnotGR$gene_biotype == "protein_coding")
  feature_annotationIdx = which(tmpAnnotGR$gene_biotype != "protein_coding")
  
  protein_coding_map = mapgro[protein_codingIdx]
  feature_annotation_map = mapgro[feature_annotationIdx]
  
  protein_codingAnnotGR = tmpAnnotGR[protein_codingIdx]
  feature_AnnotGR = tmpAnnotGR[feature_annotationIdx]
  
  protein_codingInGR = ingr[queryHits(protein_coding_map)]
  featureInGR = ingr[queryHits(feature_annotation_map)]
  
  ###############################################
  #  Calling the gene and feature annotation Functions
  ###############################################
  protCodingDF = protein_coding_AnnotationFunction_ensembl( annotGR=protein_codingAnnotGR, 
                                                            annotMap=protein_coding_map, 
                                                            inputGR=protein_codingInGR, 
                                                            inputIDs=ingr$tmpID)
  ###############################################
  #  Calling the feature annotation Function
  ###############################################
  featureDF = feature_annotation_Function_ensembl( annotGR=feature_AnnotGR, 
                                                   annotMap=feature_annotation_map, 
                                                   inputGR=featureInGR, 
                                                   inputIDs=ingr$tmpID)
  
  return(list("protCodingDF"=protCodingDF, "featureDF"=featureDF))
  
})


#' @rdname annotationSummary-method
setMethod("annotationSummary", signature("EnsemblAnnotation"), function(object, ... ){
 
 flatL = getFlatTable(object)
 
 protCodingDF = flatL$protCodingDF
 featureDF =  flatL$featureDF
  
 pcid = protCodingDF$InputID
 pcgt = protCodingDF$gene_type
 
 protCodingDFL = split(pcgt, pcid)
 
 protCodingDF_summary = protCodingDF[ which( !duplicated(protCodingDF$InputID)) ,c(1:4,6,7)]#always first entry

 protCodingDF_summary$NumberOfTranscripts =  unlist( lapply( protCodingDFL, length) )
 
 toInvestigate = protCodingDF_summary$NumberOfTranscripts > 1
 protCodingDF_summary$gene_type = as.character(protCodingDF_summary$gene_type)
 protCodingDF_summary$gene_type[toInvestigate] = unlist( lapply( protCodingDFL[toInvestigate], function(gt){
   return( paste0( na.omit(unique(gt)), collapse="/") )
 }) )
 protCodingDF_summary$gene_type[which(protCodingDF_summary$gene_type == "")] = NA
 
 rowCoverage = rowSums( featureDF[,c("feature_queryCoverage", "feature_subjectCoverage")] )# 2 is the highest number in rowSums
 rowCoverage = max(rowCoverage,na.rm=TRUE)-rowCoverage
 featureDF$rowCoverage = rowCoverage
 featureDF = featureDF[order(featureDF$InputID, featureDF$rowCoverage),]
 
 featureDF_summary = featureDF[ which( !duplicated(featureDF$InputID)),]#always first entry
 
 featid = featureDF$InputID
 featgt = featureDF$feature_type
 featureDFL = split(featgt, featid)
 
 featureDF_summary$NumberOfTranscripts =  unlist( lapply( featureDFL, function(x){ 
   return(length(x))
 }))
 
 resL = list( "protCodingDF"=protCodingDF_summary, 
              "featureDF"=featureDF_summary )
 
  return(resL)
})


