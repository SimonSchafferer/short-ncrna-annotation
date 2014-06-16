
#'Modifies the gtf file when downloaded from ensembl or ucsc by adding an 'intron' annotation line between exons
#'@param pathToGTF the path to the gtf file downloaded from UCSC or ENSEMBL
#'@param filename name of the gtf file
#'@param beginMetdata The column index indicating when the metadata starts normally column 6 (default)
#'@return RObject for annotation (UCSC, ENSEMBL) as GRanges object
#'@export
createRObject_gtf = function( pathToGTF, filename, beginMetadata=6, type="ucsc"){
  
  type = tolower(type)
  supportedTypes = c("ucsc","gencode","ensembl")
  
  if( !(type %in% supportedTypes) ){
    stop( paste0("Type not supported, this method currently supports: ", paste( supportedTypes, collapse=", ")))
  }
  
  require(rtracklayer)
  currWd = getwd()
  message(paste0("Importing file", file.path(pathToGTF,filename)) )
  gtf_file = import( file.path(pathToGTF,filename),format="gff", asRangedData=FALSE )
  gtf_file = as.data.frame( gtf_file ) 
  
  if( gencode ){
    warning("GEncode is currently not supported properly")
    #modifying some columns since files from gencode needs to have several columns separated transcript_id, gene_id separated after import by rtracklayer 
    gtfColNames = colnames(gtf_file)
    groupGTF = as.character(gtf_file)$group
    gencode_names = c("transcript_id", "gene_type", "gene_status", "gene_id", "transcript_name", "gene_name")
    transcript_type
    
    extraCols = lapply( gencode_names, function(x, groupGTF){
      res = sub(";","",sub( paste0(".*",x," " ) ,"", groupGTF))
      res = gsub("\"","",res)
      res = sub(" .*","",res)
      return(res)  
    },groupGTF )
    names(extraCols) = gencode_names
    gtf_file = gtf_file[,-which(gtfColNames == "group")]
    extraCols = do.call(cbind, extraCols)
    gtf_file = cbind(gtf_file, extraCols)
    colnames(gtf_file)[which(colnames(gtf_file) == "gene_type")] = "type" #have to change to be consistent
  } else if( ucsc ){
    #modifying some columns since files from ucsc genome browser need to have the transcript_id, gene_id separated after import by rtracklayer 
    gtf_file$transcript_id = sub(";","",sub(".*transcript_id ","",gtf_file$group))
    gtf_file$gene_id = sub(";.*","",sub("gene_id ","",gtf_file$group))
    gtf_file$transcript_id = gsub("\"","",gtf_file$transcript_id)
    gtf_file$gene_id = gsub("\"","",gtf_file$gene_id)
    gtf_file = gtf_file[,-which(colnames(gtf_file) == "group")]
  }
  
  tablename = paste0(filename,".csv")
  
  tmpOut = tempdir()
  
  write.table(gtf_file , file=file.path(tmpOut, tablename ) , sep="\t", col.names=TRUE, row.names=FALSE)
  
  modTablename =  paste0( sub(".csv","",tablename),"_mod.csv")
  tryCatch({
    setwd(file.path(tmpOut))
    system(paste0("python ",system.file("data/ModifyEnsemblAnnotation.py",package="sncRNAannotation")," -i ",
                  tablename, " -o ", modTablename))
  },error=function(e){
    setwd(currWd)
  })
   
  gtf_file = read.table( modTablename , sep="\t", header=TRUE)
  gtf_fileGR = with( gtf_file, GRanges(seqnames=seqnames, IRanges(start, end, width), strand) )
  elementMetadata(gtf_fileGR) = gtf_file[,beginMetadata:dim(gtf_file)[2]]
  
  message("Cleaning...")
  file.remove( file.path(tmpOut, tablename), file.path(tmpOut, modTablename)  )
  setwd(currWd)
  
  return( gtf_fileGR )
}

# mm9_refseq_gtf = createRObject_gtf("/home/schaffrr/test", "refSeqmm9.gtf")
# save( mm9_refseq_gtf, file="/home/schaffrr/Rworkspace/short-ncrna-annotation/data/UCSC/mm9_refseq_gtf.rda" )
# ensemblGenes = createRObject_gtf("/home/schaffrr/test", "Mus_musculus.GRCm38.75.gtf")
# hg19 = createRObject_gtf("/home/schaffrr/test", "refSeqhg19.gtf")
