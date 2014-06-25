# #Preparing libs
pathMMU = "/home/simon/dbsOfflineUse/MusMusculus/"

# refseqGenes_gtf_ucsc_mm9 = createRObject_gtf( pathMMU, "refseqGenes_gtf_ucsc_mm9.gtf", type="ucsc" )
# save(refseqGenes_gtf_ucsc_mm9, file=file.path( pathMMU, "refseqGenes_gtf_ucsc_mm9.rda" ) )
# refseqGenes_gtf_ucsc_mm10 = createRObject_gtf( pathMMU, "refseqGenes_gtf_ucsc_mm10.gtf", type="ucsc" )
# save(refseqGenes_gtf_ucsc_mm10, file=file.path( pathMMU, "refseqGenes_gtf_ucsc_mm10.rda" ) )

ensembl_gtf_v67_mm9 = createRObject_gtf( pathMMU, "ensembl_gtf_v67_mm9.gtf", type="ensembl" )
save(ensembl_gtf_v67_mm9, file=file.path( pathMMU, "ensembl_gtf_v67_mm9.rda" ) )

ensembl_gtf_v75_mm10 = createRObject_gtf( pathMMU, "ensembl_gtf_v75_mm10.gtf", type="ensembl" )
save(ensembl_gtf_v75_mm10, file=file.path( pathMMU, "ensembl_gtf_v75_mm10.rda" ) )

pathHSA = "/home/simon/dbsOfflineUse/HomoSapiens/"

# refseqGenes_gtf_ucsc_hg19 = createRObject_gtf( pathHSA, "refseqGenes_gtf_ucsc_hg19.gtf" , type="ucsc")
# save( refseqGenes_gtf_ucsc_hg19,  file=file.path( pathHSA, "refseqGenes_gtf_ucsc_hg19.rda" ))

ensembl_gtf_hg19 = createRObject_gtf( pathHSA, "ensembl_gtf_hg19.gtf", type="ensembl" )
save( ensembl_gtf_hg19,  file=file.path( pathHSA, "ensembl_gtf_hg19.rda" ))

# 
# sequence based annotation
# 
# protein_coding gene annotation
# 
# feature annotation
# 
# ensembl: features and protein_coding genes
# if no biotype is present then user must specify feature annotation or gene annotation

library(rtracklayer)
library(Biostrings)
# library(XML)
library(sncRNAannotation)

# 
# candidatesOfInterest = "/media/Rstick/workspace/ncAnnotation/data/testInterestCand.bed"
# candidatesOfInterest = import( candidatesOfInterest )
# blastcmd = "blastn -query /media/schaffrr/Rstick/workspace/ncAnnotation/data/inSeq.fsa -db /media/schaffrr/Rstick/workspace/ncAnnotation/data/testFastadb -max_target_seqs=1 -num_threads=3 -outfmt \"5\" -word_size 11"


candidateOfInterest = GRanges(seqnames=c("chr4","chr4"), IRanges( c(155429005, 155429055)  ,c(155429094,155429120)), strand=c("-","-") )
testInSeq = readDNAStringSet("/media/schaffrr/Rstick/workspace/ncAnnotation/data/inSeq.fsa")



testNew = new("NcbiBlastAnnotation", "inTest","/media/schaffrr/Rstick/workspace/ncAnnotation/data/testFastadb", testInSeq, word_size=13)

refseq = RefSeqUCSCAnnotation("refseqAnnot",system.file("resources/ucsc/", package="sncRNAannotation"),"refseqGenes_gtf_ucsc_mm9.rda", candidateOfInterest)
ensembl = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", candidateOfInterest)

ensembl = new("EnsemblAnnotation","ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", candidateOfInterest)


general = GRangesBasedAnnotation("refseqAnnot",system.file("resources/ucsc/", package="sncRNAannotation"),"refseqGenes_gtf_ucsc_mm9.rda", candidateOfInterest)

getAnnotationByID(testNew, id = 1)




qh = GRanges(seqnames="chr4", IRanges(
      start=c(1,6,1,8),
      end=c(20,14,13,19)), 
      strand="-")
sh = GRanges(seqnames="chr4", IRanges(c(5,5,5,5),c(15,15,15,15) ), strand="-")














otherAnnotationPath = "/media/Rstick/workspace/ncAnnotation/data/testAnnotation.bed"
geneAnnotPath = "/media/Rstick/workspace/ncAnnotation/data/hg19.ensembl.rda"



ensemblAnnot = new("GeneAnnotation", "ensemblAnnot",geneAnnotPath, candidatesOfInterest)
miRNAAnnot = new("GRangesBasedAnnotation", "miRNAAnnot",otherAnnotationPath, candidatesOfInterest)

#tested with blastn 2.29+ and 2.25+
# testOld = new("SequenceBasedAnnotation", "inTest","/media/Rstick/workspace/ncAnnotation/data/testFastadb", "/media/Rstick/workspace/ncAnnotation/data/inSeq.fsa")


testNew = new("NcbiBlastAnnotation", "inTest","/media/Rstick/workspace/ncAnnotation/data/testFastadb", testInSeq, word_size=13)

# convertRangesToDF(candidatesOfInterest)
# convertRangesToDF(ensemblAnnot)


bla = new("ShortRNAAnnotation", inputRanges=candidatesOfInterest[1], inputSequences=testInSeq, 
    overRanges=list("ensembl"=ensemblAnnot,"miRNA"=miRNAAnnot), overSequences=list("test"=testNew)  )

convertToUniqueAnnotationDataFrame(bla)

validObject(bla)


load("/media/Rstick/workspace/ncAnnotation/data/mm9_ensembl.rda")

#EXAMPLE protein_coding
test = mm9_ensembl[which(mm9_ensembl$gene_id == "ENSMUSG00000000702")]
testRed = reduce(test)
testSpan = GRanges( seqnames(testRed[1]), IRanges(start( testRed )[1], end(testRed)[length(testRed)]) )

intronicSeqs = gaps(testRed, start=start(testRed)[1] )
elementMetadata(intronicSeqs) = elementMetadata(test[1,])
elementMetadata(intronicSeqs)$type = "intron"
elementMetadata(intronicSeqs)$exon_number = 1:length(intronicSeqs)


test = append(test, intronicSeqs)
test = test[order(test)]

#EXAMPLE LincRNA
test = mm9_ensembl[which(mm9_ensembl$gene_id == "ENSMUSG00000085428")]
testRed = reduce(test)
testSpan = GRanges( seqnames(testRed[1]), IRanges(start( testRed )[1], end(testRed)[length(testRed)]) )

intronicSeqs = gaps(testRed, start=start(testRed)[1] )
elementMetadata(intronicSeqs) = elementMetadata(test[1,])
elementMetadata(intronicSeqs)$type = "intron"
elementMetadata(intronicSeqs)$exon_number = 1:length(intronicSeqs)
test = append(test, intronicSeqs)
test = test[order(test)]

#OK So the criteria for checking is is a GENE ID is more than once present
#EXAMPLE miRNA/snoRNA only posses one entry so no gap can be filled!
test = mm9_ensembl[which(mm9_ensembl$gene_biotype == "miRNA")][1]
testRed = reduce(test)
testSpan = GRanges( seqnames(testRed[1]), IRanges(start( testRed )[1], end(testRed)[length(testRed)]) )

intronicSeqs = gaps(testRed, start=start(testRed)[1] ) #This could result in an error please handle if gap returns an Empty GRangesObject
elementMetadata(intronicSeqs) = elementMetadata(test[1,])
elementMetadata(intronicSeqs)$type = "intron"
elementMetadata(intronicSeqs)$exon_number = 1:length(intronicSeqs)
test = append(test, intronicSeqs)
test = test[order(test)]



# ok so split the GRanges object by gene_id iterate through the list and in the end make a new list!
test = mm9_ensembl[which(mm9_ensembl$gene_id %in% c("ENSMUSG00000000702","ENSMUSG00000085428") )]

testL = split(test, test$gene_id) # results in a GRangesList
lapply( testL, function( entry ){
  
  if( length(entry) > 1 ){
    entryRed = reduce(entry)
    entrySpan = GRanges( seqnames(entryRed[1]), IRanges(start( entryRed )[1], end(entryRed)[length(entryRed)]) )
    
    intronGR = gaps(entryRed, start=start(entryRed)[1] )
    if( length(intronGR) == 0 )
    
    elementMetadata(intronicSeqs) = elementMetadata(test[1,])
    elementMetadata(intronicSeqs)$type = "intron"
    elementMetadata(intronicSeqs)$exon_number = 1:length(intronicSeqs)
    
    
  } else{
    return(entry)
  }
  


  test = append(test, intronicSeqs)
  test = test[order(test)]
  
} )


do.call( append, list(test, intronicSeqs) )


#As a second location based GRanges object: take all protein_coding genes, reduce them and fill the gaps -> naming them intergenic
#Then a region for a ncRNA is fixed: can be intergenic (between two protein_coding genes) or within a gene / intron, exon...
#All other annotations are features overlapping this region, like e.g. lincRNAs or miRNAs...







gene_biotype


#minimum: exon, CDS, start_codon, exon, CDS, exon, CDS, stop_codon, exon


hg19_v75_ensembl = import("/home/simon/dbsOfflineUse/HomoSapiens/Homo_sapiens.GRCh37.75.gtf", format="gtf", asRangedData=FALSE)

