#' @rdname xz
#' @name xzcompress
#' @title Compressing and uncompressing text files
#' 
#' @description These functions are adapted from the \code{R.utils} package from gzip to xz. Internally
#' \code{xzfile()} (see connections) is used to read (write) chunks to (from) the xz file. If the
#' process is interrupted before completed, the partially written output file is automatically removed.
#' 
#' @param filename Path name of input file.
#' @param destname Pathname of output file.
#' @param temporary If TRUE, the output file is created in a temporary directory.
#' @param skip If TRUE and the output file already exists, the output file is returned as is.
#' @param overwrite If TRUE and the output file already exists, the file is silently overwritting,
#' otherwise an exception is thrown (unless skip is TRUE).
#' @param remove If TRUE, the input file is removed afterward, otherwise not.
#' @param BFR.SIZE The number of bytes read in each chunk.
#' @param compression The compression level used (1-9).
#' @param ... Not used.
#' 
#' @return Returns the pathname of the output file. The number of bytes processed is returned as an attribute. 
#' 
#' @author Kristian Hovde Liland.
#' 
#' @examples
#' # Creating small file
#' tf <- tempfile()
#' cat(file=tf, "Hello world!")
#' 
#' # Compressing
#' tf.xz <- xzcompress(tf)
#' print(file.info(tf.xz))
#' 
#' # Uncompressing
#' tf <- xzuncompress(tf.xz)
#' print(file.info(tf))
#' file.remove(tf)
#' 
#' @export xzcompress
#' 
xzcompress <- function(filename, destname = sprintf("%s.xz", filename), temporary = FALSE, 
          skip = FALSE, overwrite = FALSE, remove = TRUE, BFR.SIZE = 1e+07, compression = 6,
          ...){
  if(!file.exists(filename)){
    stop("No such file: ", filename)
  }
  if(temporary){
    destname <- file.path(tempdir(), basename(destname))
  }
  attr(destname, "temporary") <- temporary
  if(filename == destname) 
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", filename))
  if(file.exists(destname)){
    if(skip){
      return(destname)
    } else if(!overwrite){
      stop(sprintf("File already exists: %s", destname))
    }
  }
  destpath <- dirname(destname)
  if(!file.info(destpath)$isdir) 
    dir.create(destpath)
  inn <- file(filename, open = "rb")
  on.exit(if(!is.null(inn)) close(inn))
  outComplete <- FALSE
  out <- xzfile(destname, open = "wb", compression = compression, ...)
  on.exit({
    close(out)
    if(!outComplete){
      file.remove(destname)
    }
  }, add = TRUE)
  nbytes <- 0L
  repeat{
    bfr <- readBin(inn, what = raw(0L), size = 1L, n = BFR.SIZE)
    n <- length(bfr)
    if(n == 0L) 
      break
    nbytes <- nbytes + n
    writeBin(bfr, con = out, size = 1L)
    bfr <- NULL
  }
  outComplete <- TRUE
  if(remove){
    close(inn)
    inn <- NULL
    file.remove(filename)
  }
  attr(destname, "nbrOfBytes") <- nbytes
  invisible(destname)
}
#' @rdname xz
#' @export xzuncompress
xzuncompress <- function(filename, destname = gsub("[.]xz$", "", filename, ignore.case = TRUE), 
          temporary = FALSE, skip = FALSE, overwrite = FALSE, remove = TRUE, 
          BFR.SIZE = 1e+07, ...){
  if(!file.exists(filename)){
    stop("No such file: ", filename)
  }
  if(temporary){
    destname <- file.path(tempdir(), basename(destname))
  }
  attr(destname, "temporary") <- temporary
  if(filename == destname){
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", filename))
  }
  if(file.exists(destname)){
    if(skip){
      return(destname)
    } else if(!overwrite){
      stop(sprintf("File already exists: %s", destname))
    }
  }
  destpath <- dirname(destname)
  if(!file.info(destpath)$isdir) 
    dir.create(destpath)
  inn <- xzfile(filename, open = "rb")
  on.exit(if(!is.null(inn)) close(inn))
  outComplete <- FALSE
  out <- file(destname, open = "wb")
  on.exit({
    close(out)
    if (!outComplete) {
      file.remove(destname)
    }
  }, add = TRUE)
  nbytes <- 0L
  repeat{
    bfr <- readBin(inn, what = raw(0L), size = 1L, n = BFR.SIZE)
    n <- length(bfr)
    if(n == 0L) 
      break
    nbytes <- nbytes + n
    writeBin(bfr, con = out, size = 1L)
    bfr <- NULL
  }
  outComplete <- TRUE
  if(remove){
    close(inn)
    inn <- NULL
    file.remove(filename)
  }
  attr(destname, "nbrOfBytes") <- nbytes
  invisible(destname)
}
