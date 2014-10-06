
Short ncRNA annotation package
-------------

The annotation package is written in the R programming language and. 
It provides annotation of commonly used interval/range based data such as 
the browser extended display format [(BED)](http://www.ucsc...), or the general feature format 
[(GFF)](http://www.ucsc...). In addition, it provides sequence based annotation, 
by employing the NCBI [blast+](http://www.ncbi.nlm.nih.gov/books/NBK1762/) software. 

## Interval based annotation
Basically, the genomic coordinates of the data being investigated must be given as input. This can be done by importing any of the above mentioned 
interval based file formats (GFF, BED, etc.). 

```R
library(rtracklayer)
# testInput = import( PATHTOFILE, asRangedData=FALSE)
#or simply create your own GRanges based file from a csv filee
candidateOfInterestMouse = GRanges(seqnames=c("chr4","chr4","chr7","chr5","chr16","chr4"), 
                              IRanges( c(155429005, 155429055, 20283118, 23362617,84714182,155421393),
                                      c(155429094,155429120, 20283194, 23362826, 84714202,155421564)),
                              strand=c("-","-", "+", "+","+","-") )
```

The interval based annotation can be split in two different types: 
### Feature annotation
The feature annotation can be employed for any interval based annotation. It reports all features that overlap a given range from the data that is being investigated. Therefore tracks from the UCSC genome browser may be employed such as the miRNA track, or the tRNA track. 

```R
#DO NOT EXECUTE 
general = GRangesBasedAnnotation("miRNA.annotation",system.file("resources/PATHTORESOURCE/", package="sncRNAannotation"),"mir aseBedFile.bed", candidateOfInterest)
```

### Protein-coding gene-based annotation
The protein-coding gene based annotation is based on a two-step annotation process. First the protein-coding gene is annotated separately from 
the ncRNA feature to define the reference locus. It reports the gene name strand and the overlapping position i.e. intergenic, intron, exon, 
exon/intron boundaryor UTR depending on the location of the gene that is investigated. In the second step a feature annotation is applied, 
reporting overlapping ncRNA annotations in this region i.e. the gene of interest represents snoRNA X that is located in intron Y of gene Z. 
As this annotation is more specialized, it needs to be in a format interpretable by the program. Currently the mouse and human ENSEMBL and RefSeq 
annotation is supported, respectively. However, other species may be easily employed by downloading and converting the GFF annotation from ENSEMBL, 
or the RefSeq annotation from UCSC.

Based on this data the interval based annotation searches for overlapping genomic coordinates within the annotation files provided. 
```R
library(sncRNAannotation)
ensemblMouse = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"ensembl_gtf_v67_mm9.rda", candidateOfInterestMouse)
```

A different organism may be used by downloading a gff file from [ensembl](http://...) and executing following code:
```R
#new_species_ensemblAnnot = createRObject_gtf( pathToGTF, filename, beginMetadata=6, type="ensembl")
#save(new_species_ensemblAnnot, file=PATHTOFILE/resources/ensembl/new_species_ensemblAnnot.rda)
#... and it can employed by:
#ensemblMouse = EnsemblAnnotation("ensemblAnnot", system.file("resources/ensembl/", package="sncRNAannotation"),"new_species_ensemblAnnot.rda", candidateOfI#nterestMouse)
```

The program provides a method to summarize the result to a table that features one annotation per input sample, based on the highest overlap. 
```R
testMouse = annotationSummary(ensemblMouse)
head(testMouse$protCodingDF)
head(testMouse$featureDF)
```



## Sequence based annotation
For sequence based annotation the input sequences and the path to the annotation database have to be provided. 
The program executes the NCBI blast+ with the XML output option and reports the best annotation based on the eValue. 
Currently, the precompiled frnadb database for mouse and human is providedby the package. However, other sequence database files can be employed as well. 
For conversion of fasta files to a [blast+](http://www.ncbi.nlm.nih.gov/books/NBK1762/) database files.

```R
library(Biostrings)
testseq = DNAStringSet(c( c("seq1"="CGTGTTCACAGCGGACCTTGATT"), c("seq2"="GGGCCTCAGCGTAATCCTATTGAGCATGAATTTTATGTCATCACCTTTCTTCATGA") ))
testseq_annot = NcbiBlastAnnotation("inTest", paste0(system.file("resources/frnadb/", package="sncRNAannotation"),"Mus_musculus.fasta"), testseq)
annotationSummary(testseq_annot)

```

## UML
![alt tag](https://github.com/SimonSchafferer/short-ncrna-annotation/blob/master/inst/shortNCRNAannotation.png)



