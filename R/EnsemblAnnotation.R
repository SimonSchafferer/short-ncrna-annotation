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

#' @rdname annotationSummary-method
setMethod("annotationSummary", signature("EnsemblAnnotation"), function(object, ... ){
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
 protCodingDF = protein_coding_AnnotationFunction( annotGR=protein_codingAnnotGR, 
                                    annotMap=protein_coding_map, 
                                    inputGR=protein_codingInGR, 
                                    inputIDs=ingr$tmpID)
 protCodingDFL = split(protCodingDF, protCodingDF$InputID)
 #return list and annotation summary data.frame
 protCodingDF_summary = do.call(rbind,lapply( protCodingDFL, function(x){
   df = x[1,c(1:4,6,7)]
   df$NumberOfTranscripts = dim(x)[1]
   return(df)
 }  ))
 
 ###############################################
 #  Calling the feature annotation Function
 ###############################################
 featureDF = feature_annotation_Function( annotGR=feature_AnnotGR, 
                                                   annotMap=feature_annotation_map, 
                                                   inputGR=featureInGR, 
                                                   inputIDs=ingr$tmpID)
 featureDFL = split(featureDF, featureDF$InputID)
 featureDF_summary = do.call(rbind, lapply(featureDFL ,function(x){
   #testing for completely empty subject
   test = na.omit(x)
   numberOfFeatures = dim(test)[1]
   x$NumberOfFeatures = numberOfFeatures
   if(dim(test)[1] == 0){
     return(x[1,])
   }
   
   rowCoverage = rowSums( x[,c("feature_queryCoverage", "feature_subjectCoverage")] )
   maxRowCov = which( rowCoverage == max(rowCoverage) )
   if( length(maxRowCov) > 1  ){
     return(x[1,])
   } else{
     return(x[maxRowCov,])
   }
   
 }))
 
  return( list( "protCodingDFL"=protCodingDFL, "protCodingDF"=protCodingDF_summary "featureDFL"featureDFL,"featureDF"=featureDF_summary )  )
})


protein_coding_AnnotationFunction = function( annotGR, annotMap, inputGR, inputIDs){
      
  #Subsetting to necessary information
  geneDF = with(annotGR, data.frame( 
    "InputID"=inputGR$tmpID, 
    "gene_biotype"=gene_biotype,
    "gene_type"=type,
    "gene_id"=gene_id,
    "transcript_id"=transcript_id,
    "exonIntron_number"=exon_number,
    "gene_name"=gene_name,
    "gene_transcript_name"=transcript_name
    ) )
  
  missingIDs = inputIDs[which(!inputIDs %in% geneDF$InputID)]
  
  if( length(missingIDs) != 0 ){
    missingDF = data.frame( 
      "InputID"=missingIDs, 
      "gene_biotype"=NA,
      "gene_type"="intergenic",
      "gene_id"=NA,
      "transcript_id"=NA,
      "exonIntron_number"=NA,
      "gene_name"=NA,
      "gene_transcript_name"=NA
    )
    
    geneDF = rbind( geneDF, missingDF)
  }
  rownames(geneDF) = 1:dim(geneDF)[1]
  geneDF = geneDF[order(geneDF$InputID),]
  
  return(geneDF)
}

feature_annotation_Function = function(annotGR, annotMap, inputGR, inputIDs){
  
  alignmentCoverage_featureDF = calculateAlignmentCoverageTwoGRanges( qh=inputGR, sh=annotGR )
  colnames(alignmentCoverage_featureDF) = c("feature_queryCoverage","feature_subjectCoverage")
  
  #Subsetting to necessary information
  featureDF = with(annotGR, data.frame( 
    "InputID"=inputGR$tmpID, 
    "feature_biotype"=gene_biotype,
    "feature_type"=type,
    "feature_id"=gene_id,
    "feature_transcript_id"=transcript_id,
    "feature_exonIntron_number"=exon_number,
    "feature_name"=gene_name
  ) )
  featureDF = cbind(featureDF, alignmentCoverage_featureDF)
  
  missingIDs = inputIDs[which(!inputIDs %in% featureDF$InputID)]
  
  if( length(missingIDs) != 0 ){
    missingDF = data.frame( 
      "InputID"=missingIDs, 
      "feature_biotype"=NA,
      "feature_type"=NA,
      "feature_id"=NA,
      "feature_transcript_id"=NA,
      "feature_exonIntron_number"=NA,
      "feature_name"=NA,
      "feature_queryCoverage"=NA,
      "feature_subjectCoverage"=NA
    )
    featureDF = rbind( featureDF, missingDF)
  }
  rownames(featureDF) = 1:dim(featureDF)[1]
  featureDF = featureDF[order(featureDF$InputID),]
  
  return(featureDF)
  
}



