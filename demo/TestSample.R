library(rtracklayer)
setwd("/home/simon/dbsOfflineUse/MusMusculus/")
load("/home/simon/dbsOfflineUse/MusMusculus/ensembl_gtf_v75_mm10.rda")
mm9 = import("/home/simon/dbsOfflineUse/MusMusculus/ensembl_gtf_v75_mm10.gtf",asRangedData=FALSE)
setwd("/home/simon/dbsOfflineUse/HomoSapiens")
load("/home/simon/dbsOfflineUse/HomoSapiens/ensembl_gtf_hg19.rda")
hg19 = import("/home/simon/dbsOfflineUse/HomoSapiens/ensembl_gtf_hg19.gtf", asRangedData=FALSE)
load("/home/simon/RWorkspace/short-ncrna-annotation/inst/resources/ensembl/ensembl_gtf_hg19.rda")

library(rtracklayer)
library(Biostrings)
library(sncRNAannotation)

options(stringsAsFactors=FALSE)

numberOfEntries = 1
candidateOfInterestMouse = GRanges(seqnames=rep( c("chr4","chr4","chr7","chr5","chr16","chr4"),numberOfEntries), 
                              IRanges( rep(c(155429005, 155429055, 20283118, 23362617,84714182,155421393),numberOfEntries)  ,
                                      rep(c(155429094,155429120, 20283194, 23362826, 84714202,155421564),numberOfEntries)),
                              strand=rep(c("-","-", "+", "+","+","-"),numberOfEntries) )

candidateOfInterestHum = GRanges(seqnames=rep( c("chr19","chr19","chr1"),numberOfEntries), 
                                 IRanges( rep(c(45411384, 45410998,12613 ),numberOfEntries)  ,
                                          rep(c(45411419, 45411060,12721),numberOfEntries)),
                                 strand=rep(c("+","-","+"),numberOfEntries) )

ensemblHum = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_hg19.rda", candidateOfInterestHum)
ensemblMouse = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", candidateOfInterestMouse)

generalRefSeq = GRangesBasedAnnotation("refseqAnnot",system.file("resources/ucsc/", package="sncRNAannotation"),"refseqGenes_gtf_ucsc_mm9.rda", candidateOfInterestMouse)
generalEnsembl = GRangesBasedAnnotation("test", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", candidateOfInterestMouse)
generalRefSeqHum = GRangesBasedAnnotation("testHuman", system.file("resources/ucsc/", package="sncRNAannotation"),"refseqGenes_gtf_ucsc_hg19.rda", candidateOfInterestHum)

# Rprof() 
testHum = annotationSummary(ensemblHum)
# summaryRprof() 

# Rprof() 
testMouse = annotationSummary(ensemblMouse)
# summaryRprof() 

# Rprof() 
testMgeneralRefSeq = annotationSummary(generalRefSeq)
# summaryRprof() 

# Rprof() 
testgeneralEnsembl = annotationSummary(generalEnsembl)
# summaryRprof() 
testgeneralRefSeqHum = annotationSummary(generalRefSeqHum)


library(Biostrings)
testseq = DNAStringSet(c( c("seq1"="CGTGTTCACAGCGGACCTTGATT"), c("seq2"="GGGCCTCAGCGTAATCCTATTGAGCATGAATTTTATGTCATCACCTTTCTTCATGA") ))
# testNew = new("NcbiBlastAnnotation", "inTest", paste0(system.file("resources/frnadb/", package="sncRNAannotation"),"Mus_musculus.fasta"), testseq, word_size=13)
testNew = NcbiBlastAnnotation("inTest", paste0(system.file("resources/frnadb/", package="sncRNAannotation"),"Mus_musculus.fasta"), testseq, word_size=13)
annotationSummary(testNew)
