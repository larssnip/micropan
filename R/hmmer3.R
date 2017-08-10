#' @name hmmerScan
#' @title Scanning a profile Hidden Markov Model database
#' 
#' @description Scanning FASTA formatted protein files against a database of pHMMs using the HMMER3 software.
#' 
#' @param in.files A character vector of file names.
#' @param db The full name of the database to scan.
#' @param out.folder The name of the folder to put the result files.
#' @param threads Number of CPU's to use.
#' @param verbose Logical indicating if textual output should be given to monitor the progress.
#' 
#' @details The HMMER3 software is purpose-made for handling profile Hidden Markov Models (pHMM)
#' describing patterns in biological sequences (Eddy, 2008). This function will make calls to the
#' HMMER3 software to scan FASTA files of proteins against a pHMM database. 
#' 
#' The files named in \samp{in.files} must contain FASTA formatted protein sequences. These files
#' should be prepared by \code{\link{panPrep}} to make certain each sequence, as well as the file name,
#' has a GID-tag identifying their genome. The database named in \samp{db} must be a HMMER3 formatted
#' database. It is typically the Pfam-A database, but you can also make your own HMMER3 databases, see
#' the HMMER3 documentation for help.
#' 
#' \code{\link{hmmerScan}} will query every input file against the named database. The database contains
#' profile Hidden Markov Models describing position specific sequence patterns. Each sequence in every
#' input file is scanned to see if some of the patterns can be matched to some degree. Each input file
#' results in an output file with the same GID-tag in the name. The result files give tabular output, and
#' are plain text files. See \code{\link{readHmmer}} for how to read the results into R.
#' 
#' Scanning large databases like Pfam-A takes time, usually several minutes per genome. The scan is set
#' up to use only 1 cpu per scan by default. By increasing \code{threads} you can utilize multiple CPUs, typically
#' on a computing cluster.
#' Our experience is that from a multi-core laptop it is better to start this function in default mode
#' from mutliple R-sessions. This function will not overwrite an existing result file, and multiple parallel
#' sessions can write results to the same folder.
#' 
#' @return This function produces files in the folder specified by \samp{out.folder}. Existing files are
#' never overwritten by \code{\link{hmmerScan}}, if you want to re-compute something, delete the
#' corresponding result files first.
#' 
#' @references Eddy, S.R. (2008). A Probabilistic Model of Local Sequence Alignment That Simplifies
#' Statistical Significance Estimation. PLoS Computational Biology, 4(5).
#' 
#' @note The HMMER3 software must be installed on the system for this function to work, i.e. the command
#' \samp{system("hmmscan")} must be recognized as a valid command if you run it in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panPrep}}, \code{\link{readHmmer}}.
#' 
#' @examples 
#' \dontrun{
#' # Using a FASTA file in the micropan package
#' # We need to uncompress it first...
#' extdata.path <- file.path(path.package("micropan"),"extdata")
#' filenames <- "Mpneumoniae_309_GID2.fsa"
#' pth <- lapply( file.path( extdata.path, paste( filenames, ".xz", sep="" ) ), xzuncompress )
#' 
#' # Using a miniature pHMM database in the micropan package
#' # We need to uncompress its datafiles first...
#' db <- "microfam0.hmm"
#' pth <- lapply( file.path( extdata.path,
#'    paste( db, c(".h3f.xz",".h3i.xz",".h3m.xz",".h3p.xz"), sep="" ) ), xzuncompress )
#'   
#' # ...and scanning the FASTA-file against microfam0...
#' hmmerScan(in.files=file.path(extdata.path,filenames), 
#'    db=file.path(extdata.path,db),out.folder=".")
#'   
#' # ...and compressing all files again...
#' pth <- lapply( file.path( extdata.path, filenames ), xzcompress )
#' pth <- lapply( file.path( extdata.path,
#'    paste( db, c(".h3f",".h3i",".h3m",".h3p"), sep="" ) ), xzcompress )
#' }
#' 
#' @export
hmmerScan <- function( in.files, db, out.folder, threads=0, verbose=TRUE ){
  basic <- paste( "hmmscan -o hmmer3.txt --cut_ga --noali --cpu", threads )
  db.name <- rev( unlist( strsplit( db, split="/" ) ) )[1]
  for( i in 1:length( in.files ) ){
    gi <- gregexpr( "GID[0-9]+", in.files[i], extract=T )
    rname <- paste( gi, "_vs_", db.name, ".txt", sep="" )
    res.files <- dir( out.folder )
    if( !(rname %in% res.files) ){
      if( verbose ) cat( "hmmerScan: Scanning", in.files[i], "...\n" )
      command <- paste( basic, "--domtblout", file.path( out.folder, rname ), db, in.files[i]  )
      print( command )
      system( command )
      file.remove( "hmmer3.txt" )
    }
  }
}



#' @name readHmmer
#' @title Reading results from a HMMER3 scan
#' 
#' @description Reading a text file produced by \code{\link{hmmerScan}}.
#' 
#' @param hmmer.file The name of a \code{\link{hmmerScan}} result file.
#' @param e.value Numeric threshold, hits with E-value above this are ignored (default is 1.0).
#' @param use.acc Logical indicating if accession numbers should be used to identify the hits.
#' 
#' @details The function reads a text file produced by \code{\link{hmmerScan}}. By specifying a smaller
#' \samp{e.value} you filter out poorer hits, and fewer results are returned. The option \samp{use.acc}
#' should be turned off (FALSE) if you scan against your own database where accession numbers are lacking.
#' 
#' @return The results are returned in a \samp{data.frame} with columns \samp{Query}, \samp{Hit},
#' \samp{Evalue}, \samp{Score}, \samp{Start}, \samp{Stop}, \samp{Description}. \samp{Query} is the tag
#' identifying each sequence in each genome, typically \samp{GID111_seq1}, \samp{GID121_seq3}, etc.
#' \samp{Hit} is the name or accession number for a pHMM in the database describing patterns. The
#' \samp{Evalue} is the \samp{ievalue} in the HMMER3 terminology. The \samp{Score} is the HMMER3 score for
#' the match between \samp{Query} and \samp{Hit}. The \samp{Start} and \samp{Stop} are the positions
#' within the \samp{Query} where the \samp{Hit} (pattern) starts and stops. \samp{Description} is the
#' description of the \samp{Hit}.
#' 
#' There is one line for each hit. 
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{hmmerScan}}, \code{\link{hmmerCleanOverlap}}, \code{\link{dClust}}.
#' 
#' @examples # See the example in the Help-file for dClust
#' 
#' @export
readHmmer <- function( hmmer.file, e.value=1, use.acc=TRUE ){
  al <- readLines( hmmer.file )
  al <- al[which( !grepl( "\\#", al ) )]
  al <- gsub( "[ ]+", " ", al )
  lst <- strsplit( al, split=" " )
  if( use.acc ){
    hit <- sapply( lst, function(x){ x[2] } )
  } else {
    hit <- sapply( lst, function(x){ x[1] } )
  }
  query <- sapply( lst, function(x){ x[4] } )
  ievalue <- as.numeric( sapply( lst, function(x){ x[13] } ) )
  score <- as.numeric( sapply( lst, function(x){ x[14] } ) )
  start <- as.numeric( sapply( lst, function(x){ x[18] } ) )
  stopp <- as.numeric( sapply( lst, function(x){ x[19] } ) )
  desc <- sapply( lst, function(x){ paste( x[23:length( x )], collapse=" " ) } )
  hmmer.table <- data.frame( Query=query, Hit=hit, Evalue=ievalue, Score=score, Start=start, Stop=stopp, Description=desc, stringsAsFactors=F )
  hmmer.table <- hmmer.table[which( hmmer.table$Evalue <= e.value ),]
  return( hmmer.table )
}

