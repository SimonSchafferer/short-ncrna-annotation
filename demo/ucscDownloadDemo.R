
#Saving all annotation libraries as R objects!
# gr = GRanges(seqnames="chr8",IRanges(start=9760876, end=9761003))
# gr =  GRanges(seqnames="chr1",IRanges(start=1, end=249250621))
# 
# groUCSC = fetchGeneAnnotationFromUCSC("hg19", "RefSeq",limit=5, con=con)
# 
# 
# lapply( c("ucsc","ensembl","refseq"), function(x){
#   fetchGeneAnnotationFromUCSC("hg19", x ,limit=2, con=con)
# } )

#Connect to UCSC using RMySQL
con = connectToUCSCBrowser()[["con"]]#ONLY EXECUTE ONCE!
#Wrapper function
downloadAndSaveUCSCAnnotationTables(assembly="mm9",annotationDB="refseq",path=tempdir())
#Download as GRanges object but do not save!
groUCSC = fetchGeneAnnotationFromUCSC("hg19", "refseq",limit=5, con=con)
#Download annotation as dataframe (may be quite large):
q = buildQuery( "hg19", "refGene", limit=5 )
testDF = dbGetQuery(con, q)




#Utility function library could be loaded by defining another name as the robject was named in the filesystem
#DO NOT EXECUTE
#loadLibraryAs( "anotherName", file.path(tempdir(), "hg19.refseq") )
# 
# SELECT tn.name, tn.chrom, tn.strand, tn.txStart, tn.txEnd, tn.cdsStart, tn.cdsEnd, tn.exonCount, 
# tn.exonStarts, tn.exonEnds,kginfo.category, xr.description, xr.kgID, xr.mRNA, xr.geneSymbol, xr.refseq  
# FROM mm9.refGene tn LEFT OUTER JOIN mm9.kgXref xr ON ( tn.name = xr.refseq) LEFT OUTER JOIN mm9.kgTxInfo kginfo ON ( xr.kgID = kginfo.name  ) 