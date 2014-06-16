
refSeqGenes = import( "/home/schaffrr/test/refSeqmm9.gtf",format="gff", asRangedData=FALSE )


refSeqGenes = as.data.frame( refSeqGenes ) 
refSeqGenes$transcript_id = sub(";","",sub(".*transcript_id ","",refSeqGenes$group))
refSeqGenes$gene_id = sub(";.*","",sub("gene_id ","",refSeqGenes$group))
refSeqGenes$transcript_id = gsub("\"","",refSeqGenes$transcript_id)
refSeqGenes$gene_id = gsub("\"","",refSeqGenes$gene_id)
refSeqGenes = refSeqGenes[,-which(colnames(refSeqGenes) == "group")]


write.table(refSeqGenes , file=file.path(pathToGTFfile, "refSeqGenes_gtf_mm9.csv") , sep="\t", col.names=TRUE, row.names=FALSE)
setwd(file.path(pathToGTFfile))
system(paste0("python /home/schaffrr/Rworkspace/short-ncrna-annotation/data/ModifyEnsemblAnnotation.py -i ",
              "refSeqGenes_gtf_mm9.csv", " -o ", "refSeqGenes_gtf_mm9_mod.csv"))

refSeqGenes = read.table( "refSeqGenes_gtf_mm9_mod.csv", sep="\t", header=TRUE)
refSeqGenesGR = with( refSeqGenes, GRanges(seqnames=seqnames, IRanges(start, end, width), strand) )
elementMetadata(refSeqGenesGR) = refSeqGenes[,6:dim(refSeqGenes)[2]]



#Download the ENSEMBL gtf file: 

pathToGTFfile = "/home/schaffrr/test/"
nameGTFfile = "Mus_musculus.GRCm38.75.gtf"
nameCSVfile = sub("\\.gtf","",paste0(nameGTFfile,".csv") )
library(rtracklayer)
options(stringsAsFactors=FALSE)
gtffile = import(file.path(pathToGTFfile, nameGTFfile), asRangedData=FALSE)

#skip this steps currently and load the already loaded gtf
load("/home/schaffrr/Rworkspace/short-ncrna-annotation/data/ensembl/ensembl_gtf_v67_mm9.rda")

write.table( as.data.frame( ensembl_gtf_v67_mm9 ), file=file.path(pathToGTFfile, nameCSVfile) , sep="\t", col.names=TRUE, row.names=FALSE)

#check sorting
table( unlist( sapply( unique(ensembl_gtf_v67_mm9$gene_id)[20000:20100], function(x,gtf){
  hits = which(x == gtf$gene_id)
  return( max(hits) - min(hits) + 1 == length(hits) )
}, ensembl_gtf_v67_mm9) ) )


setwd(file.path(pathToGTFfile))
system(paste0("python /home/schaffrr/Rworkspace/short-ncrna-annotation/data/ModifyEnsemblAnnotation.py -i ",
              nameCSVfile, " -o ", nameCSVfile,"_mod"))

mm9Mouse = read.table( file.path(pathToGTFfile, paste0(nameCSVfile,"_mod")), sep="\t", header=TRUE )

mm9MouseGR = with( mm9Mouse, GRanges(seqnames=seqnames, IRanges(start, end, width), strand) )
elementMetadata(mm9MouseGR) = mm9Mouse[,6:dim(mm9Mouse)[2]]
ensembl_gtf_v67_mm9_mod = mm9MouseGR
rm(mm9MouseGR)
save( ensembl_gtf_v67_mm9_mod, file="/home/schaffrr/Rworkspace/short-ncrna-annotation/data/ensembl/ensembl_gtf_v67_mm9_mod.rda" )


ensembl_gtf_v67_mm9 = ensembl_gtf_v67_mm9_mod
rm(ensembl_gtf_v67_mm9_mod)

load("/home/schaffrr/Rworkspace/short-ncrna-annotation/data/ensembl/ensembl_gtf_v67_mm9.rda")

/home/schaffrr/Rworkspace/short-ncrna-annotation/data/ensembl/ensembl_gtf_v67_mm9.rda



#minimum information needed from metadata: 

#source ( processed_pseudogene, non-coding, Mt_tRNA, etc... )
#type (exon, CDS, start_codon, gene, Selenocysteine, transcript, UTR)
#exon_number (exon number)
#gene_biotype (antisense, miRNA, pseudogene, snRNA, etc...)


library(rtracklayer)
inputTest = import("/home/schaffrr/Rworkspace/short-ncrna-annotation/data/TestData/testAnnotation.bed", format="bed")
# test = import("/home/schaffrr/gencode.v12.annotation.gtf", format="gtf", asRangedData=FALSE)
inputTest = GRanges( seqnames=c("chr4"), IRanges(155429773,155429877), name="mmu-mir-200b")

annot = subsetByOverlaps( refSeqGenesGR, inputTest )

#refSeqGenesGR
#ensembl_gtf_v67_mm9
findOverlaps( ensembl_gtf_v67_mm9[which( ensembl_gtf_v67_mm9$gene_id == "ENSMUSG00000029074" )] , inputTest)

annot = subsetByOverlaps( ensembl_gtf_v67_mm9, inputTest )

annotationMap <- findOverlaps( inputTest, annot )




#mapping the basicAnnotation file to the subset

if( length(bedToInvestigate.sub) == 0 ){
  bedToInvestigate.sub <- new("GRanges")
  annotationMap <- new("Hits")
}
