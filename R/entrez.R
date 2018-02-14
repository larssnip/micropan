#' @name entrezDownload
#' @title Downloading genome data
#' 
#' @description Retrieving genomes from NCBI using the Entrez programming utilities.
#' 
#' @param accession A character vector containing a set of valid accession numbers at the NCBI
#' Nucleotide database.
#' @param out.file Name of the file where downloaded sequences should be written in FASTA format.
#' @param verbose Logical indicating if textual output should be given during execution, to monitor
#' the download progress.
#' 
#' @details The Entrez programming utilities is a toolset for automatic download of data from the
#' NCBI databases, see \href{https://www.ncbi.nlm.nih.gov/books/NBK25500/}{E-utilities Quick Start}
#' for details. \code{\link{entrezDownload}} can be used to download genomes from the NCBI Nucleotide
#' database through these utilities.
#' 
#' The argument \samp{accession} must be a set of valid accession numbers at NCBI Nucleotide, typically
#' all accession numbers related to a genome (chromosomes, plasmids, contigs, etc). For completed genomes,
#' where the number of sequences is low, \samp{accession} is typically a single text listing all accession
#' numbers separated by commas. In the case of some draft genomes having a large number of contigs, the
#' accession numbers must be split into several comma-separated texts. The reason for this is that Entrez
#' will not accept too many queries in one chunk (less than 500). 
#' 
#' The downloaded sequences are saved in \samp{file} on your system. This will be a FASTA formatted file,
#' and should by convention have the filename extension \samp{.fsa}. Note that all downloaded sequences end
#' up in this file. If you want to download multiple genomes, you call \code{\link{entrezDownload}} multiple
#' times.
#' 
#' @return The name of the resulting FASTA file is returned (same as \code{file}), but the real result of
#' this function is the creation of the file itself.
#' 
#' @author Lars Snipen and Kristian Liland.
#' 
#' @seealso \code{\link{getAccessions}}, \code{\link[microseq]{readFasta}}.
#' 
#' @examples 
#' # Accession numbers for the chromosome and plasmid of Buchnera aphidicola, strain APS
#' acc <- "BA000003.2,AP001071.1"
#' tf <- tempfile(pattern="Buchnera_aphidicola",fileext=".fasta")
#' txt <- entrezDownload(acc,out.file=tf)
#' 
#' # Reading file to inspect
#' genome <- readFasta(tf)
#' summary(genome)
#' 
#' # ...cleaning...
#' file.remove(tf)
#' 
#' @export entrezDownload
#' 
entrezDownload <- function( accession, out.file, verbose=TRUE ){
  if( verbose ) cat( "Downloading genome..." )
  connect <- file( out.file, open="w" )
  for( j in 1:length( accession ) ){
    adr <- paste( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", accession[j], "&retmode=text&rettype=fasta", sep="" )
    entrez <- url( adr, open="rt" )
    if( isOpen( entrez ) ){
      lines <- readLines( entrez )
      writeLines( lines, con=connect )
      close( entrez )
    } else {
      cat( "Download failed: Could not open connection\n" )
    }
  }
  close( connect )
  if( verbose ) cat( "...sequences saved in", out.file, "\n" )
  return( out.file )
}


#' @name getAccessions
#' @title Collecting contig accession numbers
#' 
#' @description Retrieving the accession numbers for all contigs from a master record GenBank file.
#' 
#' @param master.record.accession The accession number (single text) to a master record GenBank file having
#' the WGS entry specifying the accession numbers to all contigs of the WGS genome.
#' 
#' @details In order to download a WGS genome (draft genome) using \code{\link{entrezDownload}} you will
#' need the accession number of every contig. This is found in the master record GenBank file, which is
#' available for every WGS genome. \code{\link{getAccessions}} will extract these from the GenBank file and
#' return them in the apropriate way to be used by \code{\link{entrezDownload}}.
#' 
#' @return A character vector where each element is a text listing the accession numbers separated by commas.
#' Each vector element will contain no more than 500 accession numbers, see \code{\link{entrezDownload}}
#' for details on this. The vector returned by \code{\link{getAccessions}} is typically used as input to
#' \code{\link{entrezDownload}}.
#' 
#' @author Lars Snipen and Kristian Liland.
#' 
#' @seealso \code{\link{entrezDownload}}.
#' 
#' @examples 
#' # The master record accession for the WGS genome Mycoplasma genitalium, strain G37
#' acc <- getAccessions("AAGX00000000")
#' # Then we use this to download all contigs and save them
#' tf <- tempfile(fileext=".fasta")
#' txt <- entrezDownload(acc,out.file=tf)
#' 
#' # Reading the file to inspect it
#' genome <- readFasta(tf)
#' summary(genome)
#' 
#' # ...cleaning...
#' file.remove(tf)
#' 
#' @importFrom microseq gregexpr
#' 
#' @export getAccessions
#' 
getAccessions <- function( master.record.accession ){
  adrId    <- paste( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=", master.record.accession, sep="" )
  idSearch <- url( adrId, open="rt" )
  if( isOpen( idSearch )){
    idDoc  <- readLines( idSearch )
    idLine <- which(regexpr("<Id>+",idDoc)==1)[1]
    id     <- substr(idDoc[idLine], 5, nchar(idDoc[idLine])-5)
    close( idSearch )
  } else {
    cat( "Download failed: Could not open connection\n" )
    return(NULL)
  }

  adr <- paste( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", id, "&retmode=text&rettype=gb", sep="" )
  entrez <- url( adr, open="rt" )
  accessions <- ""
  if( isOpen( entrez ) ){
    lines <- readLines( entrez )
    close( entrez )
    wgs.line <- gsub( "WGS[ ]+", "", lines[grep( "WGS   ", lines )] )
    ss <- unlist( strsplit( wgs.line, split="-" ) )
    head <- gregexpr( "[A-Z]+[0]+", ss[1], extract=T )
    ss.num <- as.numeric( gsub( "^[A-Z]+[0]+", "", ss ) )
    if( length( ss.num ) > 1 ){
      range <- ss.num[1]:ss.num[2]
    } else {
      range <- ss.num
    }
    ns <- ceiling( length( range ) / 500 )
    accessions <- character( ns )
    for( j in 1:ns ){
      s1 <- (j-1)*500 + 1
      s2 <- min( j*500, length( range ) )
      accessions[j] <- paste( paste( head, range[s1]:range[s2], sep="" ), collapse="," )
    }
  } else {
    cat( "Download failed: Could not open connection\n" )
  }
  return( accessions )
}
