#' @title 'Loads an R object by a given name into the globalEnvironment
#'
#' @param name Object name one wants to define
#' @param filename Path to file and name
#' @param envir environment were the object should be loaded to
#' @return TRUE
#' @docType methods
#' @export
loadLibraryGloballyAs = function( name, filename ){
  envir=globalenv()
  nameTmp = load(filename)  
  command = paste0("assign(\"",substitute(name),"\",",nameTmp,", envir = envir)")
  eval(parse(text=( command )))
  message(paste0(nameTmp, " was assigned to variable ", substitute(name)))
  return(TRUE)
}

#'@title Loads an R object by a given name and returns it
#'
#'@param name Object name one wants to define
#'@param filename Path to file and name
#'@param envir environment were the object should be loaded to
#'@return TRUE
#'@docType methods
#'@export
loadLibraryLocallyAs = function( name, filename ){
  envir=environment()
  nameTmp = load(filename)  
  command = paste0("assign(\"",substitute(name),"\",",nameTmp,", envir = envir)")
  eval(parse(text=( command )))
  return(eval(parse(text=( name ))))
}

#'@title Calculates the coverage of two GRanges objects
#'
#'@param subject GRanges object
#'@param query Granges object
#'@return data.frame coverage, once query is used as reference and once subject
#'@docType methods
#'@export
calculateAlignmentCoverageTwoGRanges = function( qh, sh ) {
  if( length(qh) != length(sh) ){
    stop("Query and Subject must be the same length!")
  }
  coverageqh = ifelse( start(qh) >= start(sh) & end(qh) <= end(sh), TRUE, FALSE ) #query is vollständig in subject enthalten
  overlapping = ifelse( start(qh) < start(sh) & end(qh) > end(sh), (end(sh)-start(sh))+1, (end(qh)-start(sh))+1 ) #subject ist vollständig in query enthalten und query steht am linken ende von subject über
  coveragesh = ifelse( start(qh) > start(sh) & end(qh) > end(sh), TRUE, FALSE ) #query steht am rechten ende von subject über
  
  overlapping[which(coverageqh)] = width(qh[which(coverageqh)])  #kompletter overlab von query in subject
  overlapping[which(coveragesh)] = end(sh[which(coveragesh)]) - start(qh[which(coveragesh)]) + 1 #rechten ende von subject über    
  
  coverageOverqh = round( overlapping/width(qh), digits = 2) #alignment coverage von query
  coverageOversh = round( overlapping/width(sh), digits = 2)#alignment coverage von subject  
  #return the alignment coverage of each Granges object to each object always taking the longer covering one (relative to length) as reference
  return(data.frame("queryCoverage"=coverageOverqh, "subjectCoverage"=coverageOversh))
}


