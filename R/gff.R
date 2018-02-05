#' @name gff2fasta
#' @title Retrieving sequences from genome
#' 
#' @description Retrieving the sequences specified in a \code{gff.table}.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' @param genome A \code{\link{Fasta}} object with the genome sequence(s).
#' 
#' @details Each row in \code{gff.table} (see \code{\link{readGFF}}) describes a genomic feature
#' in the \code{genome}. The information in the columns Seqid, Start, End and Strand are used to retrieve
#' the sequences from \code{genome$Sequence}. Every Seqid in the \code{gff.table}
#' must match the first token in one of the \code{genome$Header} texts.
#' 
#' @return A \code{\link{Fasta}} object with one row for each row in \code{gff.table}. 
#' The \code{Header} for each sequence is a summary of the information in the
#' corresponding row of \code{gff.table}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{findOrfs}}.
#' 
#' @examples Mer her
#' 
#' @useDynLib micropan
#' @importFrom Rcpp evalCpp
#' 
#' @export gff2fasta
#' 
gff2fasta <- function( gff.table, genome ){
  utag <- unique( gff.table$Seqid )
  tagz <- unique( sapply( strsplit( genome$Header, split=" " ), function(x){x[1]} ) )
  if( sum( is.na( match( utag, tagz ) ) ) > 0 )
    stop( "Seqid's in the gff.table do not match the genome$Header" )
  
  n <- nrow( gff.table )
  fobj <- data.frame( Header=gffSignature( gff.table ), 
                      Sequence=rep( "", n ),
                      stringsAsFactors=F )
  for( i in 1:length( tagz ) ){
    idx <- which( gff.table$Seqid == tagz[i] )
    seq <- extractSeq( genome$Sequence[i], 
                       gff.table$Start[idx], 
                       gff.table$End[idx], 
                       gff.table$Strand[idx] )
    idd <- which( gff.table$Strand[idx] < 0 )
    seq[idd] <- reverseComplement( seq[idd] )
    fobj$Sequence[idx] <- seq
  }
  class( fobj ) <- c( "Fasta", "data.frame" )
  return( fobj )
}



#' @name gffSignature
#' @title GFF signature text
#' 
#' @description Making a signature text from \code{gff.table} data.
#' 
#' @param gff.tableat A \code{gff.table} (\code{data.frame}) with genomic features information.
#' 
#' @details For each row in \code{gff.table} a text is created by pasting these
#' data together, adding some explanatory text. This function is used by \code{link{gff2fasta}}
#' to create the Header-lines for the \code{\link{Fasta}} object when retrieving the sequences.
#' 
#' @return A vector of texts, one for each row in \code{gff.table}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{findOrfs}}, \code{\link{gff2fasta}}.
#' 
#' @examples Mer her
#' 
#' @useDynLib micropan
#' 
#' @export gffSignature
#' 
gffSignature <- function( gff.table ){
  desc <- paste( "Seqid=", gff.table$Seqid,
                 ";Source=", gff.table$Source,
                 ";Type=", gff.table$Type,
                 ";Start=", gff.table$Start,
                 ";End=", gff.table$End,
                 ";Score=", gff.table$Score,
                 ";Strand=", gff.table$Strand,
                 ";Phase=", gff.table$Phase,
                 ";Attributes=", gff.table$Attributes,
                 sep="" )
  return( desc )
}



#' @name readGFF
#' @aliases readGFF writeGFF
#' @title Reading and writing GFF-tables
#' 
#' @description Reading or writing a \code{gff.table} from/to file.
#' 
#' @param in.file Name of file with a GFF-table.
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' @param out.file Name of file.
#' 
#' @details A \code{gff.table} is simply a \code{data.frame} with columns
#' adhering to the format specified by the GFF3 format, see
#' https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details. There is
#' one row for each feature.
#' 
#' The following columns should always be in a full \code{gff.table} of the GFF3 format:
#' \itemize{
#'   \item Seqid. A unique identifier of the genomic sequence on which the feature resides.
#'   \item Source. A description of the procedure that generated the feature, e.g. \code{"R-package micropan::findOrfs"}.
#'   \item Type The type of feature, e.g. \code{"ORF"}, \code{"16S"} etc.
#'   \item Start. The leftmost coordinate. This is the start if the feature is on the Sense strand, but
#'   the end if it is on the Antisense strand.
#'   \item End. The rightmost coordinate. This is the end if the feature is on the Sense strand, but
#'   the start if it is on the Antisense strand.
#'   \item Score. A numeric score (E-value, P-value) from the \code{Source}. 
#'   \item Strand. A \code{"+"} indicates Sense strand, a \code{"-"} Antisense.
#'   \item Phase. Only relevant for coding genes. the values 0, 1 or 2 indicates the reading frame, i.e. 
#'   the number of bases to offset the \code{Start} in order to be in the reading frame.
#'   \item Attributes. A single string with semicolon-separated tokens prociding additional information.
#' }
#' Missing values are described by \code{"."} in the GFF3 format. This is also done here, except for the
#' numerical columns Start, End, Score and Phase. Here \code{NA} is used, but this is replaced by
#' \code{"."} when writing to file. 
#' 
#' @return \code{readGFF} returns a \code{gff.table} with the columns described above.
#' 
#' \code{writeGFF} writes the supplied \code{gff.table} to a text-file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfTable}}, \code{\link{lorf}}.
#' 
#' @examples Mer her
#' 
#' @export readGFF writeGFF
#' 
readGFF <- function( in.file ){
  gff.table <- read.table( in.file, sep="\t", header=F,
                           stringsAsFactors=F, comment.char="#" )
  if( ncol( gff.table ) != 9 ) stop( "File", in.file, "does not contain data in GFF3 format" )
  colnames( gff.table ) <- c( "Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes" )
  w <- options()$warn
  options( warn=-1 )
  gff.table$Start <- as.numeric( gff.table$Start )
  gff.table$End <- as.numeric( gff.table$End )
  gff.table$Score <- as.numeric( gff.table$Score )
  gff.table$Phase <- as.numeric( gff.table$Phase )
  options(warn=w)
  return( gff.table )
}
writeGFF <- function( gff.table, out.file ){
  line1 <- c("##gff-version 3.2.1")
  lines <- sapply( 1:nrow(gff.table), function(i){paste( gff.table[i,], collapse="\t" )} )
  lines <- gsub( "\tNA\t", "\t.\t", lines )
  lines <- gsub( "\tNA$", "\t.", lines )
  writeLines( c( line1, lines ), con=out.file )
  return( NULL )
}
