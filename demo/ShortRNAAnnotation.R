# TODO: Add comment
# 
# Author: simon
###############################################################################
#
#library("XML")
#library("rtracklayer")
#library("Biostrings")

##################################################################################################################
#
#	Extended Chip Annotation class (supports csv, bed, gtf and ncbi blast xml files) 
#
##################################################################################################################

setClass("ShortRNAAnnotation", 
representation( inputRanges="GRanges", inputSequences="XStringSet", overRanges="list", overSequences="list"  )) 

validShortRNAAnnotation=function(object) {
  
  if( length(object@inputRanges) != length(object@inputSequences) ){
    return("Length of InputRanges must match the length of the inputSequences")
  }
  
    allRanged = unlist( lapply( object@overRanges, function(x){
      return(is(x,"RangedAnnotation"))
    } ) )
    
  if( sum(!allRanged) != 0 ){
      return("Only RangedAnnotation objects are allowed in slot: OverRanges")
  }

  allSequence = unlist( lapply( object@overSequences, function(x){
      return(is(x,"SequenceBasedAnnotation"))
    } ) )
    if( sum(!allRanged) != 0 ){
      return("Only SequenceBasedAnnotation objects are allowed in slot: OverSequences")
    }
  return(TRUE)
}

setValidity("ShortRNAAnnotation", validShortRNAAnnotation)

#initialize will need a GRanges object, a DNAStringSet and a ID

######################
#	Annotation GRanges and DNAString Set should be in the same order and have the same size and have an internal ID column
######################
setMethod("initialize", "ShortRNAAnnotation", function(.Object, inputRanges,inputSequences, overRanges, overSequences){
      .Object@inputRanges = inputRanges
      .Object@inputSequences = inputSequences
			.Object@overRanges = overRanges
			.Object@overSequences = overSequences
			return(.Object)
		})

ShortRNAAnnotation = function( inputRanges,inputSequences, overRanges, overSequences ){
  new(  "ShortRNAAnnotation", inputRanges,inputSequences, overRanges, overSequences )
}

setMethod("show", "ShortRNAAnnotation", function(object){	
#   ira = slot(object, "inputRanges")
#   print( length(ira)  )
#   
#   iseq = slot(object, "inputSequences")
#   print( length(iseq) )
#   
#   overR = slot(object, "overRanges")
#   print( length(overR) )
#   
#   overSeq = slot(object, "overSequences")
#   print( "" length(overSeq) )
  
  print(slot(object, "inputRanges"))
  print(slot(object, "inputSequences"))
  print(slot(object, "overRanges"))
  print(slot(object, "overSequences"))
  
		})

######################
#	Accessor Methods
######################
setGeneric("inputRanges", function(object) { standardGeneric("inputRanges") })
setMethod("inputRanges",signature(object="ShortRNAAnnotation"),function(object) {
			slot(object, "inputRanges")
		})

setGeneric("inputSequences", function(object) { standardGeneric("inputSequences") })
setMethod("inputSequences",signature(object="ShortRNAAnnotation"),function(object) {
			slot(object, "inputSequences")
		})
setGeneric("overRanges", function(object) { standardGeneric("overRanges") })
setMethod("overRanges",signature(object="ShortRNAAnnotation"),function(object) {
			slot(object, "overRanges")
		})
setGeneric("overSequences", function(object) { standardGeneric("overSequences") })
setMethod("overSequences",signature(object="ShortRNAAnnotation"),function(object) {
			slot(object, "overSequences")
		})
























setGeneric("printAnnotationByID", function(object, id) { standardGeneric("printAnnotationByID") })
setMethod("printAnnotationByID",signature(object="ShortRNAAnnotation"),function(object, id=NULL) {
			
			if(length(Annotations(object)) > 1 & is.null(id)){
				stop("Please provide only one candidate/ID")
			} 
			
			if(!is.null(id) & length(id) == 1){ object <- getAnnotationByID(object,id); } else{ stop("Please provide only one candidate/ID") }
			
			annotShow <- lapply(Annotations(object), function(x){ 
						if( is(x,"RangedAnnotation") ){
							
							if( !(length(annotationGR(x)) == 0 & length(annotationGRAS(x)) == 0) ){
								return(x)
							}							
						} else if( is(x,"SequenceBasedAnnotation")){
							
							if( length(BlastOutput(x)) != 0 ){
								return(x)
							}
						}
					})
			
			annotShow <- annotShow[unlist(lapply(annotShow, function(x){!is.null(x)}))]
			object@Annotations <- annotShow
			return(object)	
			
		})


#' Accessor method annotName
#' 
#' @export
#' @docType methods
#' @rdname annotName-method
setGeneric("annotName", function(object) { standardGeneric("annotName") })

######################
#	Other Methods
######################
setGeneric( "getAnnotationsAsGRanges", function(object, ... ){ standardGeneric("getAnnotationsAsGRanges") } )
setMethod( "getAnnotationsAsGRanges" ,signature(object="ShortRNAAnnotation"),function(object) {
			
			annotations <- Annotations(object)
			annotationsGR <- lapply(annotations, function(x){
						#calling getAnnotationByID with 
						
						if( is(x,"RangedAnnotation") ){
							x <- as( x ,"RangedAnnotation" )
						}
						
						gro <- annotationGR(x)
						gro <- append(gro, annotationGRAS(x) )
						
						return(gro)
						
					})
			names(annotationsGR) <- sub("_bed*","",names(Annotations(analysisContainerL$snpAnnotation)) )
			
			#reduce all to only name and score
			reduceToNameScore <- function(gro){
				
				score = elementMetadata(gro)$score
				name = elementMetadata(gro)$name
				
				if(is.null(score)){
					score = 0
				} else{
					score = as.numeric(score)
				}
				
				if( is.null( name ) ){
					#first try type
					if(!is.null(elementMetadata(gro)$type)){
						name = elementMetadata(gro)$type
					} else{
						name = "NA"
					}
				} else{
					name = as.character(name)
				}
				
				elementMetadata(gro) = NULL
				elementMetadata(gro)$score = score
				elementMetadata(gro)$name = name
				
				return(gro)
				
			}
			
			annotationsGR = lapply(annotationsGR, reduceToNameScore)
			
			variabledef = paste( "\"",names(annotationsGR),"\"","=annotationsGR[[\"",names(annotationsGR),"\"]]", sep="" )
			variabledef = paste(variabledef, collapse=",")
			
			eval(parse(text=paste("tmpvarstorage = GRangesList(",variabledef,")", sep="") ))
			
			return(tmpvarstorage)
			
		})



setGeneric( "getAnnotationByID", function(object, id, ... ){ standardGeneric("getAnnotationByID") } )
setMethod( "getAnnotationByID" ,signature(object="ShortRNAAnnotation"),function(object, id) {
			
			if(length(id) > 1){
				stop("Please provide only one ID")
			}
			
			ids <- IDMap( object )
			idIndex <- as.numeric(which( id == ids ))
			
			annotations <- Annotations(object)
			
			annotations <- lapply(annotations, function(x){
						#calling getAnnotationByID with 
						
						if( is(x,"RangedAnnotation") ){
							x <- as( x ,"RangedAnnotation" )
						}
						
						getAnnotationByID(x, idIndex)
						
					})
			
			object@Annotations <- annotations
			
			return(object)
			
		})

setGeneric("convertToDataFrame", function(object, ...){ standardGeneric("convertToDataFrame") })
setMethod("convertToDataFrame", "ShortRNAAnnotation", function(object, ...){
			
			idsToTraverse <- IDMap(object)		
			
			######################
			#	Conversion thoughts: 
			#	BED/CSV Files: 
			#		dbName_Orientation:sense/as/none ; dbName_Name:GeneName
			#	GTF Files: 
			#		dbName_Orientation:sense/as/none ; db_Source:mm9_snp128 ; db_Type:exon/intron ; groupt:gene_id
			#	Blast Files: (taking the first entry -> since prioritized by eValue)
			#		dbName_Hit_def:miRNA128 ; dbName_Hsp_align_len:22 ; dbName_Hsp_eValue:0.000045
			######################
			
			annotations <- lapply(Annotations(object), function(x){
						#calling getAnnotationByID with 
						if( is(x,"RangedAnnotation") ){
							
							x <- as( x ,"RangedAnnotation" )
							#LOOP Through all
							rangedDF <- lapply( 1:length(idsToTraverse), function(idIndex){
										
										rangedAnnot <- getAnnotationByID(x, idIndex)
										convertToDataFrame(rangedAnnot)
										
										
									} )
							
							rangedDF <- do.call(rbind,rangedDF)
							
							return(rangedDF)
							
						} else if( is(x,"SequenceBasedAnnotation") ){
							sequenceDF <-lapply( 1:length(idsToTraverse), function(idIndex){
										
										blastout <- getAnnotationByID(x, idIndex)
										convertToDataFrame( blastout )
										
									} )
							
							sequenceDF <- do.call(rbind,sequenceDF)
							return(sequenceDF)
						}
						
					})
			names(annotations) <- NULL 
			return( do.call(cbind, annotations) )
			
		})


setGeneric("convertToUniqueAnnotationDataFrame", function(object, ...){ standardGeneric("convertToUniqueAnnotationDataFrame") })
setMethod("convertToUniqueAnnotationDataFrame", "ShortRNAAnnotation", function(object, ...){
			
			idsToTraverse <- IDMap(object)		
			
			######################
			#	Output Table should contain
			#	ncRNAType	ncRNAName	ncRNAStrandness nearestGene	positionInGene	strandnessOfGene	SeqIdentity
			#
			#
			#	Conversion thoughts: 
			#	BED/CSV Files: 
			#		dbName_Orientation:sense/as/none ; dbName_Name:GeneName
			#	GTF Files: 
			#		dbName_Orientation:sense/as/none ; db_Source:mm9_snp128 ; db_Type:exon/intron ; groupt:gene_id
			#	Blast Files: (taking the first entry -> since prioritized by eValue)
			#		dbName_Hit_def:miRNA128 ; dbName_Hsp_align_len:22 ; dbName_Hsp_eValue:0.000045
			######################
			
			annotations <- lapply(Annotations(object), function(x){
						#calling getAnnotationByID with 
						if( is(x,"RangedAnnotation") ){
							
							x <- as( x ,"RangedAnnotation" )
							#LOOP Through all
							rangedDF <- lapply( 1:length(idsToTraverse), function(idIndex){
										
										rangedAnnot <- getAnnotationByID(x, idIndex)
										convertToDataFrame(rangedAnnot)
									} )
							
							rangedDF <- do.call(rbind,rangedDF)
							
							return(rangedDF)
							
						} else if( is(x,"SequenceBasedAnnotation") ){
							sequenceDF <-lapply( 1:length(idsToTraverse), function(idIndex){
										blastout <- getAnnotationByID(x, idIndex)
										convertToDataFrame( blastout )
									} )
							sequenceDF <- do.call(rbind,sequenceDF)
							return(sequenceDF)
						}
						
					})
			names(annotations) <- NULL 
			return( do.call(cbind, annotations) )
			
		})

setGeneric("ncRNAAnnotation", function(object, groupdf, ...){ standardGeneric("ncRNAAnnotation") })
setMethod("ncRNAAnnotation", "ShortRNAAnnotation", function(object, groupdf=data.frame(), groupNames=c(), ...){
			
			idsToTraverse <- IDMap(object)		
			#groups should be based on the annotations that should be coerced -> possible in bed files by coverage, in blastFiles by e-value in gene files by priority (e.g. RefSeq -> ENSEMBL -> UCSC)
			#therefore this is dataframe with a priority and group column -> based on this the annotation is ordered and sorted!
			groupdfCols = c("Annotation", "Groups")
			if( !sum( colnames(groupdf) %in% groupdfCols ) == length(groupdfCols) ){
				warning("Please Check the group data.frame columns, should be Annotation and Groups")
				stop()
			}
			
			if( !length(groupNames) == length(unique(groupdf$Groups)) ) {
				warning("Groups will be named automatically");
				groupNames = unique(groupdf$Groups)
			}
			
			
			if(dim( groupdf ) [1] == 0) {
				groupdf = data.frame("Annotation"=c(1:length(Annotations(object))) , "Groups" = c(1:length(Annotations(object))) )
			}
			
			######################
			#	Conversion thoughts: 
			#	BED/CSV Files: 
			#		dbName_Orientation:sense/as/none ; dbName_Name:GeneName
			#	GTF Files: 
			#		dbName_Orientation:sense/as/none ; db_Source:mm9_snp128 ; db_Type:exon/intron ; groupt:gene_id
			#	Blast Files: (taking the first entry -> since prioritized by eValue)
			#		dbName_Hit_def:miRNA128 ; dbName_Hsp_align_len:22 ; dbName_Hsp_eValue:0.000045
			######################			
			annotations = lapply( unique(groupdf$Groups), function( g ) {
						
						annotIndex = groupdf[groupdf$Groups %in% g,]$Annotation
						
						#Check if valid groups are formed: Always Ranged Annotations and Blast annotations in one group! 
						if(!( sum( unlist( lapply( Annotations(object)[annotIndex], inherits, "RangedAnnotation" ) ) ) == length(annotIndex) | sum( unlist( lapply( Annotations(object)[annotIndex], inherits, "RangedAnnotation" ) ) ) == 0 )){
							warning("The groups that have been defined are not from the same class -> always group bed with bed and blast with blast! Currently only blast gene and bed grouping is allowed")
							stop()
						} 
						
						#merging the groups in the different classes
						annotationsToMerge = Annotations(object)[annotIndex]
						
						if( is(annotationsToMerge[[1]],"RangedAnnotation") ){
							
							######################
							#	internal method calculates the coverage in percent of two GRanges objects
							######################
							calculateAlignmentCoverageTwoGRangesGenes = function( sh, qh ) {
								
								#calculating the amount of coverage in percent (both from subject and query)
								coverageqh = ifelse( start(qh) >= start(sh) & end(qh) <= end(sh), TRUE, FALSE ) #query is vollständig in subject enthalten
								overlapping = ifelse( start(qh) < start(sh) & end(qh) > end(sh), (end(sh)-start(sh))+1, (end(qh)-start(sh))+1 ) #subject ist vollständig in query enthalten und query steht am linken ende von subject über
								coveragesh = ifelse( start(qh) > start(sh) & end(qh) > end(sh), TRUE, FALSE ) #query steht am rechten ende von subject über
								overlapping[coverageqh] = width(qh[coverageqh])  #kompletter overlab von query in subject
								overlapping[coveragesh] = end(sh[coveragesh]) - start(qh[coveragesh]) + 1 #rechten ende von subject über
								
								coverageOverqh = round( overlapping/width(qh), digits = 2) #alignment coverage von query
								coverageOversh = round( overlapping/width(sh), digits = 2)#alignment coverage von subject	
								
								if( "blocks" %in% colnames(elementMetadata(sh)) ){
									ncRNAVect = unlist( lapply( elementMetadata(sh)$blocks, function(x){ length(x) == 1 } ) )
								} else{
									ncRNAVect = rep(TRUE, length(sh))
								}
								
								#in case of refseq  modify the gene selection on exons
								refseqncRNAs = grep( "[NR,XR]_.*", elementMetadata(sh)$name )
								if( length(refseqncRNAs) != 0 ) {
									ncRNAVect[refseqncRNAs] = TRUE #in case of refseq
								}
								refseqProteinMRNA = grep( "[NM,NP,XP,XM]_.*", elementMetadata(sh)$name )
								if( length(refseqProteinMRNA) != 0 ) {
									ncRNAVect[refseqProteinMRNA] = FALSE #in case of refseq
								}
								#in case of ensembl modify the gene selection on exons
								ensemblncRNA = grep( "ENS.{,3}T.*", elementMetadata(sh)$name )#transcript
								if( length(ensemblncRNA) != 0 ) {
									ncRNAVect[ensemblncRNA] = TRUE #in case of ensembl
								}
								ensemblProt = grep( "ENS.{,3}[E,G,P].*", elementMetadata(sh)$name )
								if( length(ensemblProt) != 0 ) {
									ncRNAVect[ensemblProt] = FALSE #in case of ensembl
								}
								
								blocks = elementMetadata(sh)$blocks #the width and end is mixed up
								
								#This will take some time unfortunately!!!
								positionInGene = mapply( function( queryIndex, block ){
#									print(queryIndex)				
#									print(length(block))
											
											tmpIRanges = IRanges( start = width(block), end = start(block) )
											qhCurrent = qh[queryIndex]
											
											overlapsExons = findOverlaps(ranges(qhCurrent), tmpIRanges)
											
											if( length(overlapsExons) > 0 ){
												#overlap in exon found!
												#Coverage from query
												subjectHit = tmpIRanges[subjectHits(overlapsExons)]
												
												coverage = ifelse( start(qhCurrent) >= start(subjectHit) & end(qhCurrent) <= end(subjectHit), TRUE, FALSE ) #query is vollständig in subject enthalten
												
												if( sum(coverage) > 0 ){ #it could be that the RNA spans over two exons with a small intron -> then coverage contains more than one value and it would be exon/intron 
													return( "exon" )
												}else{
													return("exon/intron")
												}
											} else{
												return("intron")
											}
											
										}, as.list(c(1:length(blocks))), as.list(blocks) )
								
								
								
								#return the alignment coverage of each Granges object to each object always taking the longer covering one (relative to length) as reference
								return(data.frame("qh"=coverageOverqh, "sh"=coverageOversh, "ncRNA"=ncRNAVect, "GeneRegion"=positionInGene))
							}
							
							calculateAlignmentCoverageTwoGRanges = function( sh, qh ) {
								coverageqh = ifelse( start(qh) >= start(sh) & end(qh) <= end(sh), TRUE, FALSE ) #query is vollständig in subject enthalten
								overlapping = ifelse( start(qh) < start(sh) & end(qh) > end(sh), (end(sh)-start(sh))+1, (end(qh)-start(sh))+1 ) #subject ist vollständig in query enthalten und query steht am linken ende von subject über
								coveragesh = ifelse( start(qh) > start(sh) & end(qh) > end(sh), TRUE, FALSE ) #query steht am rechten ende von subject über
								overlapping[coverageqh] = width(qh[coverageqh])  #kompletter overlab von query in subject
								overlapping[coveragesh] = end(sh[coveragesh]) - start(qh[coveragesh]) + 1 #rechten ende von subject über
								
								coverageOverqh = round( overlapping/width(qh), digits = 2) #alignment coverage von query
								coverageOversh = round( overlapping/width(sh), digits = 2)#alignment coverage von subject	
								
								if( "blocks" %in% colnames(elementMetadata(sh)) ){
									ncRNAVect = unlist( lapply( elementMetadata(sh)$blocks, function(x){ length(x) == 1 } ) )
								} else{
									ncRNAVect = rep(TRUE, length(sh))
								}
								
								#in case of refseq  modify the gene selection on exons
								refseqncRNAs = grep( "[NR,XR]_.*", elementMetadata(sh)$name )
								if( length(refseqncRNAs) != 0 ) {
									ncRNAVect[refseqncRNAs] = TRUE #in case of refseq
								}
								refseqProteinMRNA = grep( "[NM,NP,XP,XM]_.*", elementMetadata(sh)$name )
								if( length(refseqProteinMRNA) != 0 ) {
									ncRNAVect[refseqProteinMRNA] = FALSE #in case of refseq
								}
								#in case of ensembl  modify the gene selection on exons
								ensemblncRNA = grep( "ENS.{,3}T.*", elementMetadata(sh)$name )#transcript
								if( length(ensemblncRNA) != 0 ) {
									ncRNAVect[ensemblncRNA] = TRUE #in case of ensembl
								}
								ensemblProt = grep( "ENS.{,3}[E,G,P].*", elementMetadata(sh)$name )
								if( length(ensemblProt) != 0 ) {
									ncRNAVect[ensemblProt] = FALSE #in case of ensembl
								}
								
								#return the alignment coverage of each Granges object to each object always taking the longer covering one (relative to length) as reference
								return(data.frame("qh"=coverageOverqh, "sh"=coverageOversh, "ncRNA"=ncRNAVect))
							}
							
							
							######################
							#	Special Case of Gene annotation -> but is always further processed in general range based annotation
							# 	subset annotationsToMerge in this case
							######################
							annotationsToMergeGene = annotationsToMerge[ unlist( lapply( annotationsToMerge, is, "GENEDFAnnotation" ) ) ]
							bestAnnotDFGene = data.frame()
							prefix <- Name((annotationsToMerge[[1]]) )
							
							if( length(annotationsToMergeGene) > 0 ){
								
								##########################################################
								#	Gene Annotation!! non coding RNAs are marked as: txEnd == cdsStart == cdsEnd -> no exon intron structure! (longes = gene)
								#	Output: ID (candidate), GeneID, geneName, Region, Strandness, Distance, ncRNA 
								#	
								#	REFSEQ
								#	refseqID = new.env()
								#	refseqID[["Complete_genomic_molecule"]] = "NC_.*"
								#	refseqID[["Incomplete_genomic_region"]] = "NG_.*"
								#	refseqID[["mRNA"]] = "NM_.*"
								#	refseqID[["ncRNA"]] = "NR_.*"
								#	refseqID[["Protein"]] = "NP_.*"
								#	refseqID[["predicted_mRNA_model"]] = "XM_.*"
								#	refseqID[["predicted_ncRNA_model"]] = "XR_.*"
								#	refseqID[["predicted_Protein_model"]] = "XP_.*"
								#	
								#	ENSEMBL
								#	ensemblID = new.env()
								#	ensemblID[["exon"]] = ".{1,6}E.*"
								#	ensemblID[["protein_family"]] = ".{1,6}FM.*"
								#	ensemblID[["gene"]] = ".{1,6}G.*"
								#	ensemblID[["gene_tree"]] = ".{1,6}GT.*"
								#	ensemblID[["protein"]] = ".{1,6}P.*"
								#	ensemblID[["regulatory_feature"]] = ".{1,6}R.*"
								#	ensemblID[["transcript"]] = ".{1,6}T.*"
								##########################################################
								
								mapOfBedAnnotationGene = lapply( annotationsToMergeGene, function( group ) {
																					
											qh = InputRanges(object)[ queryHits( annotationMap( group ) ) ]
											sh = annotationGR(group)[ subjectHits( annotationMap( group ) ) ]
											
											if(length(qh) == 0){
												#return an empty data.frame
												return( data.frame( "queryHits"=1:length(InputRanges(object)), 
																"subjectHits"=rep(NA,length(InputRanges(object))), 
																"qh"=rep(NA, length(InputRanges(object))), 
																"sh"=rep(NA, length(InputRanges(object))), 
																"ncRNA"=rep(FALSE, length(InputRanges(object))),
																"GeneRegion"=rep(NA, length(InputRanges(object))),
																"meanAlignCoverage"=rep(NA, length(InputRanges(object))),
																"GeneLength"=rep(NA, length(InputRanges(object)))
														)
												)
											}
											
											coverageOverDF = calculateAlignmentCoverageTwoGRangesGenes( "sh"=sh, "qh"=qh )
											#get the coverage level of this annotation -> in order to compare with other annotations in this group and decide which annotation has more coverage and therefore is better suited
											coverageOverDF$meanAlignCoverage = rowMeans(coverageOverDF[,c("sh","qh")] )
											coverageOverDF$GeneLength = width(sh)
											
											#if the alignement coverage is on the oposite strand -> subtract 1 from the result then it will no change the order if all annotations are on the oposite strand, but
											#if a better coverage is on the other strand it will be scored as worse than on the same strand!
											coverageOverDF$meanAlignCoverage = ifelse( as.character(strand(sh)) == as.character(strand(qh)) | as.character(strand(qh)) == "*", coverageOverDF$meanAlignCoverage, coverageOverDF$meanAlignCoverage-1 )
											
											
											map = as.data.frame(annotationMap( group )) 
											map = cbind(map, coverageOverDF)
											
											indexBestCov = sapply( unique(map$queryHits), function(x){
														index = which( map$queryHits == x )
														index = index[ which( map[index,]$GeneLength == max(map[index,]$GeneLength) & !map[index,]$ncRNA ) ]
														return( index[1] ) 
													} )
											
											
											mapTmp = data.frame("queryHits"=1:length(InputRanges(object)))
											map = merge( mapTmp, map[indexBestCov,], by="queryHits", sort=FALSE, all.x=TRUE )
											
											map = map[order(map$queryHits), ]
											
											rownames(map) = 1:dim(map)[1]
											
											return( map[order(map$queryHits),] )
										})
								inputRanges = lapply( annotationsToMergeGene, function( group ) {
											return( InputRanges(object)[ queryHits( annotationMap( group ) ) ] )
										})
								
								evals = lapply( mapOfBedAnnotationGene, function( x ) {
											return(x$meanAlignCoverage)
										})
								evals = as.data.frame(do.call( cbind, evals ))
								
								chooselistToTraverse = apply( evals, 1, function( x ) { 
											if( sum(is.na(x)) == length(x) ){
												return(NA) 
											} else{ 
												return(which( !is.na(x))[1]) 
											} 
										}  )
								traverseL = lapply( c(1: length(annotationsToMergeGene)), function(x){ return( which( x == chooselistToTraverse) ) }) 
								
								annot = annotationsToMergeGene[[1]];trav = traverseL[[1]]; moba=mapOfBedAnnotationGene[[1]];inputRange= inputRanges[[1]]
								
								bestAnnotL = mapply(  function(annot, trav, moba, inputRange){
											if( length(trav) != 0 ){
												index = as.integer(moba[trav,"subjectHits"])
#											print(index)
												currentAnnotationGR = annotationGR(annot)[index]
												currentInput = inputRange[index]
												name = elementMetadata( currentAnnotationGR )[,1] 
												strand = as.character(strand(currentAnnotationGR))
												#intron or intergenic, distance to gene if NA Gene Name if possible, length of gene
												#ok Find out about UTR coding Sequence exon intron ...
												df = data.frame( "ID" = trav, "GeneName" = name, "GeneAlignCov"=moba[trav,"meanAlignCoverage"], 
														"GeneStrand" = strand, "GeneRegion"=moba[trav,"GeneRegion"], 
														"GeneLength"=moba[trav,"GeneLength"], "AnnotRef"=rep(Name(annot), length(trav) ) ) 
												return( df )									
											}
										},annotationsToMergeGene,  traverseL, mapOfBedAnnotationGene, inputRanges, SIMPLIFY=FALSE )
								
								bestAnnotDF = do.call(rbind, bestAnnotL)
								
								tmpdf = data.frame("ID" = rownames(evals))
								
								bestAnnotDFGene = merge( bestAnnotDF, tmpdf, by="ID", all.y=TRUE, sort=FALSE)
								bestAnnotDFGene = bestAnnotDFGene[order(as.numeric(bestAnnotDFGene$ID)),]
								colnames(bestAnnotDFGene)[1] = "GeneID"
								rownames(bestAnnotDFGene) = 1:dim(bestAnnotDFGene)[1]
								#modify the annotationsToMerge that is left!
							} 
							
							######################
							#	General ranged annotation -> take the annotation that has the highest coverage! 
							######################							
							mapOfBedAnnotation = lapply( annotationsToMerge, function( group ) {
										
										
										qh = InputRanges(object)[ queryHits( annotationMap( group ) ) ]
										sh = annotationGR(group)[ subjectHits( annotationMap( group ) ) ]
										
										if(length(qh) == 0){
											#return an empty data.frame
											emptyDF = data.frame( "queryHits"=1:length(InputRanges(object)), 
													"subjectHits"=rep(NA,length(InputRanges(object))), 
													"qh"=rep(NA, length(InputRanges(object))), 
													"sh"=rep(NA, length(InputRanges(object))), 
													"ncRNA"=rep(FALSE, length(InputRanges(object))),
													"GeneRegion"=rep(NA, length(InputRanges(object))),
													"meanAlignCoverage"=rep(NA, length(InputRanges(object))),
													"GeneLength"=rep(NA, length(InputRanges(object)))
											)
											return( emptyDF)
										}
										
										coverageOverDF = calculateAlignmentCoverageTwoGRanges( sh, qh )
										#get the coverage level of this annotation -> in order to compare with other annotations in this group and decide which annotation has more coverage and therefore is better suited
										coverageOverDF$meanAlignCoverage = rowMeans(coverageOverDF[,c("sh","qh")] )
										
										#if the alignement coverage is on the oposite strand -> subtract 1 from the result then it will no change the order if all annotations are on the oposite strand, but
										#if a better coverage is on the other strand it will be scored as worse than on the same strand!
										coverageOverDF$meanAlignCoverage = ifelse( as.character(strand(sh)) == as.character(strand(qh)) | as.character(strand(qh)) == "*", coverageOverDF$meanAlignCoverage, coverageOverDF$meanAlignCoverage-1 )
										
										map = as.data.frame(annotationMap( group )) 
#								map$meanAlignCoverage = coverageOverDF$meanAlignCoverage
										map = cbind(map, coverageOverDF)
										
										
										indexBestCov = sapply( unique(map$queryHits), function(x){
													index = which( map$queryHits == x )
													index = index[ which( map[index,]$meanAlignCoverage == max(map[index,]$meanAlignCoverage)  & map[index,]$ncRNA  ) ]
													return( index[1] ) 
												} )
										
										indexBestCov = na.omit(indexBestCov)
										if( length( indexBestCov ) == 0){
											emptyDF = data.frame( "queryHits"=1:length(InputRanges(object)), 
													"subjectHits"=rep(NA,length(InputRanges(object))), 
													"qh"=rep(NA, length(InputRanges(object))), 
													"sh"=rep(NA, length(InputRanges(object))), 
													"ncRNA"=rep(FALSE, length(InputRanges(object))),
													"GeneRegion"=rep(NA, length(InputRanges(object))),
													"meanAlignCoverage"=rep(NA, length(InputRanges(object))),
													"GeneLength"=rep(NA, length(InputRanges(object)))
											)
											return(emptyDF)
										} 
										
										
										mapTmp = data.frame("queryHits"=1:length(InputRanges(object)))
										map = merge( mapTmp, map[indexBestCov,], by="queryHits", sort=FALSE, all.x=TRUE )
										map = map[order(map$queryHits), ]
#								rownames(map) = 1:dim(map)[1]
										
										return( map[order(map$queryHits),] )
									})
							
							evals = lapply( mapOfBedAnnotation, function( x ) {
										return(x$meanAlignCoverage)
									})
							evals = as.data.frame(do.call( cbind, evals ))
							
							chooselistToTraverse = apply( evals, 1, function( x ) { if( sum(is.na(x)) == length(x) ){ return(NA) } else{ which( x == max(x, na.rm=TRUE) )[1] } }  )
							traverseL = lapply( c(1: length(annotationsToMerge)), function(x){ return( which( x == chooselistToTraverse) ) }) 
							
							bestAnnotL = mapply(  function(annot, trav, moba){
										if( length(trav) != 0 ){
											index = as.integer(moba[trav,"subjectHits"])
											name = elementMetadata( annotationGR(annot)[index] )[,1] 
											strand = as.character(strand(annotationGR(annot)[index]))
											return( data.frame( "ID" = trav, "ncRNAName" = name, "ncRNAAlignCov"=moba[trav,"meanAlignCoverage"], "ncRNAStrand" = strand, "ncAnnotRef"=rep(Name(annot), length(trav) )  ) ) #Go through mapping file!!! annotationGR(annot)[ trav ] )									
										}
									},annotationsToMerge,  traverseL, mapOfBedAnnotation, SIMPLIFY=FALSE )
							
							bestAnnotDF = do.call(rbind, bestAnnotL)
							
							if( is.null(bestAnnotDF) ){	
								#NOTHING FOUND IN ALL THE ANNOTATIONS!
								bestAnnotDF = data.frame( "ncID" = rownames(evals), "ncRNAName" = rep(NA,length(rownames(evals))), 
										"ncRNAAlignCov"=rep(NA,length(rownames(evals))), "ncRNAStrand" = rep(NA,length(rownames(evals))), 
										"ncAnnotRef"=rep(NA,length(rownames(evals))) ) #Go through mapping file!!! annotationGR(annot)[ trav ] )
								if( dim(bestAnnotDFGene)[1] != 0 ){
									bestAnnotDF = cbind(bestAnnotDFGene, bestAnnotDF[,-1])#except the ID
								}
								
								return( bestAnnotDF )
							}
							
							tmpdf = data.frame("ID" = rownames(evals))
							
							bestAnnotDF = merge( bestAnnotDF, tmpdf, by="ID", all.y=TRUE, , sort=FALSE)
							colnames(bestAnnotDF)[1] = "ncID"
							bestAnnotDF = bestAnnotDF[order(as.numeric(bestAnnotDF$ncID)),]
							rownames(bestAnnotDF) = 1:dim(bestAnnotDF)[1]
							
							if( dim(bestAnnotDFGene)[1] != 0 ){
								bestAnnotDF = cbind(bestAnnotDFGene, bestAnnotDF[,-1])#except the ID
							}
							
							return(bestAnnotDF)					
							
						}
						######################
						#	Merging Sequence based annotations
						######################
						if( is(annotationsToMerge[[1]],"SequenceBasedAnnotation") ){
							
							#take the name of the first entry as prefix > evaluate the highest e-number... and then the priority -> take first one if e-values are the same!
							
							if( length( annotationsToMerge ) > 1 ){
								evals = as.data.frame( lapply( annotationsToMerge, function( group ) {
													minEval = unlist( lapply(BlastOutput(group), function(x){ if( is.null(x$Hsp_evalue) ){t = 10}else{t = min( x$Hsp_evalue )}; return(t) } ) )
													return( as.numeric( minEval ) )
												}) )
								
								chooselistToTraverse = apply( evals, 1, function( x ) { which( x == min(x) )[1] }  )
								
								traverseL = lapply( 1:length(annotIndex), function(x){ return( which( x == chooselistToTraverse) ) }) 
								
								bestAnnotL = mapply(  function(annot, trav){
											if( length(trav) != 0 ){
												return( BlastOutput(annot)[ trav ] )									
											} 
										},annotationsToMerge,  traverseL )
								
								if( class(bestAnnotL[[1]]) == "list" | is.null(bestAnnotL[[1]]) ){
									tmpL = list()
									bestAnnotL = lapply( bestAnnotL, function( x ) { tmpL <<- append(tmpL, x) } )
									bestAnnotL = tmpL
									rm(tmpL)
								}
								names(bestAnnotL) = sub(".*Query_", "", names(bestAnnotL))
								bestAnnotL = bestAnnotL[order(as.integer(names(bestAnnotL)))]
								names(bestAnnotL) = paste("Query", names(bestAnnotL), sep="_")
								
								if( ! length(bestAnnotL) == length(BlastOutput(annotationsToMerge[[1]])) ) {
									warning("Sequence Based Annotation unification Error -> Please Check")
								}
							} else{
								bestAnnotL = BlastOutput(annotationsToMerge[[1]])
							}
							
							prefix <- DBName((annotationsToMerge[[1]]) )
							
							sequenceDF <- do.call(rbind,lapply( bestAnnotL, function( lEntry ){
												tmp <- c( lEntry[1,"Hit_def"], 
														( as.numeric(lEntry[1,"Hsp_align-len"]) - as.numeric(lEntry[1,"Hsp_gaps"]) ) / as.numeric(lEntry[1,"Hit_len"]),
														as.character( (as.numeric( lEntry[1,"Hsp_align-len"] ) - as.numeric(lEntry[1,"Hsp_gaps"]) ) / ( as.numeric(lEntry[1,"Hsp_query-to"]) - as.numeric(lEntry[1,"Hsp_query-from"]) + 1 ) ), 
														(lEntry[1,"Hsp_evalue"]) ) 
												
												if(length( tmp ) == 0 || as.numeric(tmp[4]) >  0.01) { tmp = c(NA, NA, NA, NA)} else{ tmp[2] = round(as.numeric(tmp[2]), 2); tmp[3] = round(as.numeric(tmp[3]), 2) }
												names(tmp) = c("hitdef", "hitAlign","queryAlign", "evalue" ) 
												return(tmp)
											} ) )
							
#					colnames(sequenceDF) = paste(prefix, colnames(sequenceDF), sep=".")
							sequenceDF = as.data.frame(sequenceDF)
							return(sequenceDF)
						}
						
					}) 
			
			
			names(annotations) <- groupNames 
			#do.call(cbind, annotations)
			return( annotations )
		})		

##########################################################
#	Other Methods for post-annotation
##########################################################	
annotIDenv = new.env()
annotIDenv[["ensembl_transcript_id"]] = "ENS.{,3}T.*"
annotIDenv[["ensembl_gene_id"]] = "ENS.{,3}G.*"
annotIDenv[["refseq_mrna"]] = "NM_.*"
annotIDenv[["refseq_ncrna"]] = "NR_.*"
annotIDenv[["ucsc"]] = "uc[0-9]{2,}.*"
annotIDenv[["frnadb"]] = "FR[0-9]{3,}.*"
annotIDenv[["mirbase"]] = "MI...[0-9]+"

setGeneric("extendBlastDataframe", function(df, ...){ standardGeneric("extendBlastDataframe") })
setMethod("extendBlastDataframe", signature(df="data.frame"), function(df = data.frame(), ...){
			require("RCurl")
			blastCol = "hitdef"
			#First extend the frnadb entries 
			occIndex = grep(annotIDenv[["frnadb"]], df[,blastCol] )
			frnaDBids = gsub( "\\|.*","", df[,blastCol] )
			frnaDBDescription = gsub( ".*\\|","", df[,blastCol] )
			
			frnaDBtmp = data.frame("ID"=frnaDBids, "Organism"=rep("",length(frnaDBids)), "Description"=frnaDBDescription)
			
			occIndex = grep(annotIDenv[["frnadb"]], frnaDBids )
			
			infos = sapply( as.character(frnaDBtmp[occIndex,"ID"]), fetchFRNADBOrgansism )
			
			tmp = mapply( function(info, index){
						#may be not the best approach
						frnaDBtmp[index,"Organism"] <<- info
						return(NULL)
					},infos, occIndex )
			
			#Next extend miBase entries
			occIndexMiRNA = grep(annotIDenv[["mirbase"]], df[,blastCol] )
			mirBaseids = sub( " .*","",sub( ".*-miR-.{1,4} ","", df[,blastCol] , ignore.case=TRUE) )
			mirBaseOrganism = sub( " miR.*","", sub( paste(".*",annotIDenv[["mirbase"]], " ", sep=""),"", df[,blastCol] ), ignore.case=TRUE )
			mirBaseDescription = gsub( " .*","", df[,blastCol] )
			
			frnaDBtmp$ID[occIndexMiRNA] = mirBaseids[occIndexMiRNA]
			frnaDBtmp$Organism[occIndexMiRNA] = mirBaseOrganism[occIndexMiRNA]
			frnaDBtmp$Description[occIndexMiRNA] = mirBaseDescription[occIndexMiRNA]
			
			df$hitdef = frnaDBtmp$ID
			df$Organism = frnaDBtmp$Organism
			df$Description = frnaDBtmp$Description
			
			addFRNADBBiotype = function( df ){
				df$Biotype = c("")
				grc = function( x ){ grep(x, df$Description, ignore.case=TRUE) }
				
				df$Biotype[ grc( "SRP" ) ] = "SRP"
				df$Biotype[ grc( "snoRNA" ) ] = "snoRNA"
				df$Biotype[ grc( "miR-" ) ] = "miRNA"
				df$Biotype[ grc( "tRNA" ) ] = "tRNA"
				df$Biotype[ grc( "rRNA" ) ] = "rRNA"
				df$Biotype[ grc( "scRNA" ) ] = "scRNA"
				df$Biotype[ grc( "Vault" ) ] = "Vault RNA"
				df$Biotype[ grc( "Y RNA" ) ] = "Y RNA"
				df$Biotype[ grc( "4.5S" ) ] = "4.5S RNA"
				df$Biotype[ grc( "7SK" ) ] = "7SK RNA"
				df$Biotype[ grc( "piRNA" ) ] = "piRNA"
				df$Biotype[ grc( "putative conserved" ) ] = "ncRNA prediction"
				df$Biotype[ grc( "SSU" ) ] = "SSU RNA"
				df$Biotype[ grc( "Small subunit ribosomal" ) ] = "SSU RNA"
				df$Biotype[grc("spliceosomal")] = "splicosomal RNA"
				df$Biotype[ grc( "LSU" ) ] = "LSU RNA"
				df$Biotype[ grc( "5S" ) ] = "5S RNA"
				df$Biotype[grc("scaRNA")] = "scaRNA"
				df$Biotype[grc("snRNA")] = "snRNA"
				df$Biotype[grc("snmRNA")] = "snmRNA"
				df$Biotype[ grc( "non-protein coding.*transcript")] = "processed_transcript"	
				df$Biotype = ifelse(df$Biotype == "" & !is.na(df$Description), "other", df$Biotype)
				return(df)
			}
			
			df = addFRNADBBiotype(df)
			
			return(df)
			
		}) 


setGeneric("fetchFRNADBOrgansism", function(frnadbID, ...){ standardGeneric("fetchFRNADBOrgansism") })
setMethod("fetchFRNADBOrgansism", signature(frnadbID="character"), function(frnadbID = c(), ...){
			require("RCurl")
			require("XML")
			empty = NA
			#getting through REST service to FRNADB entries then parsing the XML file: 
			xmlannot = try( xmlTreeParse(sub("<result xmlns.*<query>", "<result>\n<query>", getURL(paste("http://www.ncrna.org/frnadb/doc/api/entry/", frnadbID, sep=""))), asText = TRUE, useInternal = TRUE) )
			if( class(xmlannot)[1] == "try-error" ){
				return( empty )
			}
			
			top = xmlRoot(xmlannot)
			
			if(xmlValue(top[[1]]) == "Not Found"){
				return(empty)
			}
			
#		description = xmlValue(xmlannot[["//description"]])
			organismNames = c("Mus musculus","Homo sapiens", "Rattus norvegicus")
			#determining the organism
			
			organisms = sapply( organismNames, function(orgn){
						if( length( getNodeSet(top,paste("//result/entry_list/entry/organism/ncbi_taxonomy[@name=\'",orgn,"\']" ,sep="") ) ) != 0 ){
							return( orgn )
						} else {
							return(NA)
						}
					} )
			
			organism = if(length( organisms[which(!is.na(organisms))] ) == 0 ){ organism = "other" } else{ organisms[which(!is.na(organisms))][1] }
#		if( length( getNodeSet(top,"//result/entry_list/entry/organism/ncbi_taxonomy[@name='Mus musculus']") ) != 0 ){
#			organism = "Mus musculus"
#		} else if( length( getNodeSet(top,"//result/entry_list/entry/organism/ncbi_taxonomy[@name='Homo sapiens']") ) != 0){
#			organism = "Homo sapiens"
#		} else{
#			organism = "other"
#		}		
			return( organism )
			
		})


setGeneric("extendGeneDataframe", function(df, ...){ standardGeneric("extendGeneDataframe") })
setMethod("extendGeneDataframe", signature(df="data.frame"), function(df = data.frame(), dataset="mmusculus_gene_ensembl", ...){
			
			#ENSEMBLE Features: http://www.gencodegenes.org/gencode_biotypes.html
			#IG_C_gene, IG_D_gene, IG_J_gene, IG_V_gene, TR_C_gene, TR_J_gene, TR_V_gene, TR_D_gene, Immunoglobulin genes
			#IG_C_pseudogene, IG_J_pseudogene, IG_V_pseudogene, TR_V_pseudogene, TR_J_pseudogene, Inactivated Immunoglobulin genes
			#Mt_rRNA, Mt_tRNA, miRNA, misc_RNA, rRNA, snRNA, snoRNA -> RFAM/mirBase non-coding RNA sequences
			# Mt_tRNA_pseudogene, tRNA_pseudogene, snoRNA_pseudogene, snRNA_pseudogene, scRNA_pseudogene, rRNA_pseudogene, misc_RNA_pseudogene, miRNA_pseudogene ncRNAs predicted to be pseudogenes (by ENSEMBL)
			#TEC, nonsense_mediated_decay, non_stop_decay, retained_intron, protein_coding, processed_transcript, non_coding, ambiguous_orf, sense_intronic, sense_overlapping, antisense, pseudogene, 
			#processed_pseudogene, polymorphic_pseudogene, retrotransposed, transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene, unitary_pseudogene, unprocessed_pseudogene
			#artifact, lincRNA, LRG_gene, 3prime_overlapping_ncrna, disrupted_domain
			
			library("biomaRt")
			ensembl=useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org/biomart/martservice/", dataset=dataset)
			
			callBiomart = function( filter, ids, col ){
				
				
				bmout = getBM(attributes=c("transcript_biotype","external_transcript_id"), filters=filter,values=ids, mart=ensembl)
				index = 1
				if(dim(bmout)[1] > 1 && bmout$transcript_biotype[1] == "lincRNA"){
					bmout$transcript_biotype[1] = paste( bmout$transcript_biotype, collapse=";" )
					bmout$external_transcript_id[1] = paste( bmout$external_transcript_id, collapse=";" )						
				} else if( dim(bmout)[1] > 1 ){
					
					isncRNA=bmout$transcript_biotype %in% c("Mt_rRNA", "Mt_tRNA", "miRNA", "misc_RNA", "rRNA", "snRNA", "snoRNA")
					if(sum(isncRNA) > 0){
						index = which(isncRNA)[1]
					}
				}
				
				return( bmout[index,] )	
			}
			
			fillMartDF = function( martDF, entryID, col, df ){
				
				occIndex = grep(annotIDenv[[entryID]], df[,col] )
				ids = df[occIndex,col]
				if(length(ids) != 0){
					
					#This can be change to callBiomart( entryID, ids )  -> way faster!
					
					biomRet = lapply( ids, function(id){
								callBiomart( entryID, id, col ) 
							} )
					biomL = mapply( function(biom, occ){
								if(dim(biom)[1] != 0){
									martDF[occ,] <<- biom[1,]										
								} else{
									return(c("",""))
								}
							}, biomRet, occIndex)
					
					
				}
				return(martDF)
			}
			
			geneNameCol = "GeneName"
			ncNameCol = "ncRNAName"
			
			annotIDenv[["ensembl_transcript_id"]] = "ENS.{,3}T.*"
			annotIDenv[["ensembl_gene_id"]] = "ENS.{,3}G.*"
			annotIDenv[["refseq_mrna"]] = "NM_.*"
			annotIDenv[["refseq_ncrna"]] = "NR_.*"
			annotIDenv[["ucsc"]] = "uc[0-9]{2,}.*"
			
			######################
			#	Check for Gene Biotype
			######################
			geneDF = data.frame("GeneBiotype"=rep("",dim(df)[1]), "TranscriptID"=rep("",dim(df)[1]), stringsAsFactors=FALSE )
			
			geneDF = fillMartDF( geneDF, "refseq_mrna", geneNameCol, df )
			geneDF = fillMartDF( geneDF, "ensembl_transcript_id", geneNameCol, df )
			geneDF = fillMartDF( geneDF, "ensembl_gene_id", geneNameCol, df )
			geneDF = fillMartDF( geneDF, "refseq_mrna", geneNameCol, df )
			geneDF = fillMartDF( geneDF, "refseq_ncrna", geneNameCol, df )
			geneDF = fillMartDF( geneDF, "ucsc",geneNameCol, df )
			
			geneMod = df[,1:7]
			geneMod = cbind(geneMod, geneDF)
			
			ncRNADF = data.frame("ncRNABiotype"=rep("",dim(df)[1]), "ncRNATranscriptID"=rep("",dim(df)[1]), stringsAsFactors=FALSE )
			
			ncRNADF = fillMartDF( ncRNADF, "refseq_mrna", ncNameCol, df )
			ncRNADF = fillMartDF( ncRNADF, "ensembl_transcript_id", ncNameCol, df )
			ncRNADF = fillMartDF( ncRNADF, "ensembl_gene_id", ncNameCol, df )
			ncRNADF = fillMartDF( ncRNADF, "refseq_mrna", ncNameCol, df )
			ncRNADF = fillMartDF( ncRNADF, "refseq_ncrna",ncNameCol, df )
			ncRNADF = fillMartDF( ncRNADF, "ucsc", ncNameCol, df )
			
			ncMod = df[,8:11]
			ncMod = cbind(ncMod, ncRNADF)
			
			geneMod = cbind(geneMod, ncMod)
			
			######################
			#	lincRNA one exon correction!
			######################
			indicesDuplicate = grep( ";", geneMod$ncRNABiotype )
			#two cases -> only lincRNAs found
			#lincRNA / snoRNA found
			if( length(indicesDuplicate) != 0){
				allLincRNAs = as.logical( sapply(geneMod[indicesDuplicate,"ncRNABiotype"],function(bt){ sum(!( unlist(strsplit( bt,";")) %in% "lincRNA" ) ) == 0 }))
				geneMod[indicesDuplicate,"ncRNABiotype"] = ifelse( allLincRNAs , "lincRNA", gsub( "lincRNA;","",geneMod[indicesDuplicate,"ncRNABiotype"]) )
				
			}
			return(geneMod)
			
		}) 

setGeneric("modifyRepeatDataframe", function(df, ...){ standardGeneric("modifyRepeatDataframe") })
setMethod("modifyRepeatDataframe", signature(df="data.frame"), function(df = data.frame(), ...){
			
			if( "ncRNAName" %in% colnames(df)  ){
				dfTmp = df
				dfTmp$ncRNAFull = dfTmp$ncRNAName
				dfTmp$ncRNAName = sub("\\|.*","",df$ncRNAName)
				return(dfTmp)
			} else{
				warning("No ncRNA Name column present")
				return(df)
			}
			
		})
















# require("XML")
# require("rtracklayer")
# require("Biostrings")
# 
# if(is.null(annotGRanges)){
#   stop("GRanges Class has to be provided")
# } else if( class(annotGRanges) != "GRanges" ){
#   stop("GRanges Class has to be provided")
# } else if( is.na(genome(annotGRanges)[1]) ){
#   stop("GRangesForUCSCGenome Class has to be provided -> since the genome is required!")
# }
# 
# if(is.null(id)){
#   id <- elementMetadata(annotGRanges)$ID
# }
# if( length(id) != length(annotGRanges) ){
#   stop("Eather in the GRanges object should be an ID in the elementMetadata or an additional vector of IDs (same length as GRanges object) should be provided")
# }
# 
# #variable definitions (static)
# tmpFastaPath <- "/tmp/blastTmpFile.fasta"
# annotationFileHeader <- c("AnnotationFileName", "AnnotationType", "PathToAnnotation")
# supportedRanges <- c("bed","gtf", "csv", "genedf")
# supportedSeq <- c("blast")
# 
# #annotation model has to be provided
# if( !is.null(annotationModel) && class(annotationModel) == "character" ){
#   annotationModel <- read.table(annotationModel, sep="\t", header=TRUE, stringsAsFactors = FALSE)				
# } else if( !is.null(annotationModel) && class(annotationModel) == "data.frame" ){
#   
#   if( dim(annotationModel)[2] == length(annotationFileHeader) && sum( colnames(annotationModel) == annotationFileHeader) == length(annotationFileHeader) ){
#     annotationModel <- annotationModel
#   } else{
#     stop("Please Check your annotation data.frame! The columns, should be: ", paste(annotationFileHeader,collapse=" ") )
#   }
# } else{
#   stop("AnnotationModel is eather Null or not in supported format (character, data.frame)")
# }
# 
# 
# if( !is.null( annotDNAString ) ){
#   writeXStringSet( annotDNAString, tmpFastaPath , format="fasta" ) 
# } else if( "blast" %in% annotationModel[,annotationFileHeader[2]] ){
#   stop("Blast Search is defined in the annotation Model (csv) but no DNAStringSet was given!")
# }
# 
# 
# #Check the data types provided
# if( sum( annotationModel[,annotationFileHeader[2]] %in% c(supportedRanges, supportedSeq) ) != dim(annotationModel)[1] ){
#   stop(paste("Annotation Type not known!","\nSupported Types: ",paste(c(supportedRanges, supportedSeq),collapse=",") ) )				
# }
# 
# if( sum(colnames(annotationModel) != annotationFileHeader) != 0 ){
#   stop("Please Check your annotation file header, should be: ", paste(annotationFileHeader,collapse=" ") )
# }
# 
# .Object@Annotations <- apply( annotationModel, 1, function(ca){
#   
#   if( ca[annotationFileHeader[2]] %in% supportedSeq ){
#     #Creating a sequence based annotation object
#     print(paste("Calling Sequence based annotation for", ca[annotationFileHeader[1]] ))
#     new( "SequenceBasedAnnotation", ca[annotationFileHeader[1]], ca[annotationFileHeader[3]], tmpFastaPath )
#     
#   } else if( ca[annotationFileHeader[2]] %in% supportedRanges ){
#     #Creating a range based Annotation Object 
#     className <- paste(toupper(ca[annotationFileHeader[2]]),"Annotation",sep="")
#     print(paste("Calling Range based annotation for", ca[annotationFileHeader[1]], "as", className))
#     new( className, ca[annotationFileHeader[1]], ca[annotationFileHeader[3]], annotGRanges )
#   } else{
#     stop(paste("Annotation Type not known!","\nSupported Types: ",paste(c(supportedRanges, supportedSeq),collapse=",") ) )
#   }
#   
# } )
#names(.Object@Annotations) <- apply(annotationModel[,c(annotationFileHeader[1],annotationFileHeader[2])], 1, paste, collapse = "_")
