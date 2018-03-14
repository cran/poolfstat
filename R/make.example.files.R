#' Create example files
#' @description Write in the current directory example files corresponding to a sync (as obtained when parsing mpileup files with PoPoolation) and vcf (as obtained when parsing mpileup files with VarScan) gzipped files
#' @param writing.dir Directory where to copy example files (e.g., set writing.dir=getwd() to copy in the current working directory)
#' @examples
#'  make.example.files(writing.dir=tempdir())
#' @export
make.example.files <- function(writing.dir=""){
  if(writing.dir==""){stop("ERROR: Please provide the directory path where to copy the example files (e.g., set writing.dir=getwd() to copy in the current working directory)")}
  file.copy(from=system.file('ex.sync.gz',package='poolfstat'),to=paste0(writing.dir,'/ex.sync.gz'))
  file.copy(from=system.file('ex.vcf.gz',package='poolfstat'),to=paste0(writing.dir,'/ex.vcf.gz'))
}


