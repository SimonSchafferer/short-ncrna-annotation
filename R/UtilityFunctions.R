#'Loads an R object by a given name into the globalEnvironment
#'@param name Object name one wants to define
#'@param filename Path to file and name
#'@param envir environment were the object should be loaded to
#'@return TRUE
#'@export
loadLibraryGloballyAs = function( name, filename ){
  envir=globalenv()
  nameTmp = load(filename)  
  command = paste0("assign(\"",substitute(name),"\",",nameTmp,", envir = envir)")
  eval(parse(text=( command )))
  message(paste0(nameTmp, " was assigned to variable ", substitute(name)))
  return(TRUE)
}
#'Loads an R object by a given name and returns it
#'@param name Object name one wants to define
#'@param filename Path to file and name
#'@param envir environment were the object should be loaded to
#'@return TRUE
#'@export
loadLibraryLocallyAs = function( name, filename ){
  envir=environment()
  nameTmp = load(filename)  
  command = paste0("assign(\"",substitute(name),"\",",nameTmp,", envir = envir)")
  eval(parse(text=( command )))
  return(eval(parse(text=( name ))))
}
