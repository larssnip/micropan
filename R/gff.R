#' @name gff2fasta
#' @title Retrieving sequences from genome
#' 
#' @description Retrieving the sequences specified in a \code{gff.table}.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' @param genome A \code{\link{Fasta}}-object with the genome sequence(s).
#' 
#' @details Each row in \code{gff.table} (see \code{\link{findOrfs}}) describes a genomic feature in the \code{genome}.
#' This includes the \code{Seqid} indicating the genomic sequence, the coordinates \code{Start}
#' and \code{Stop} as well as the \code{Strand}. Every \code{Seqid} in the \code{gff.table}
#' must match the first token in one of the \code{genome$Header} texts.
#' 
#' @return A \code{\link{Fasta}} object with one row for each row in \code{gff.table}. 
#' The \code{Header} for each sequence is a summary of the \code{gff.table} information in the
#' corresponding row of \code{gff.table}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{findOrfs}}, \code{\link{lorfs}}.
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
                       gff.table$Stop[idx], 
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
                 ";Type=", gff.table$Type,
                 ";Start=", gff.table$Start,
                 ";Right=", gff.table$Stop,
                 ";Strand=", gff.table$Strand,
                 ";Attributes=", gff.table$Attributes,
                 sep="" )
  return( desc )
}



#' @name readOrfTable
#' @aliases readOrfTable writeOrfTable
#' @title Reading and writing ORF tables
#' 
#' @description Reading or writing a data.frame with ORF information from/to file.
#' 
#' @param in.file Name of file with an ORF-table.
#' @param orf.table A data.frame listing ORF information, see \code{\link{orfTable}}.
#' @param out.file Name of file.
#' 
#' @details The \code{orf.table} must be a \code{data.frame} with the columns specified in
#' \code{\link{orfTable}}. These functions use \code{\link{read.table}} and \code{\link{write.table}}
#' to read/write such a table from/to a textfile.
#' 
#' @return \code{readOrfTable} returns a data.frame with ORF information, see \code{\link{orfTable}}.
#' 
#' \code{writeOrfTable} writes an ORF-table to a text-file.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfTable}}, \code{\link{lorf}}.
#' 
#' @examples Mer her
#' 
#' @export readOrfTable writeOrfTable
#' 
readOrfTable <- function( in.file ){
  orf.table <- read.table( in.file, sep="\t", header=T,
                           colClasses=c( "character", "numeric", "numeric", "numeric", "numeric" ),
                           stringsAsFactors=F )
  return( orf.table )
}
writeOrfTable <- function( orf.table, out.file ){
  write.table( orf.table, file=out.file, sep="\t", row.names=F )
  return( NULL )
}
