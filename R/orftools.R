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



#' @name orfLength
#' @title ORF lengths
#' 
#' @description Computes the lengths of all ORFs in an ORF-table.
#' 
#' @param orf.table A data.frame listing ORF information, see \code{\link{orfTable}}.
#' 
#' @return A vector of ORF lengths, measured as the number of amino acids after translation.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfTable}}.
#' 
#' @examples Mer her
#' 
#' @export orfLength
#' 
orfLength <- function( orf.table ){
  return( (abs( orf.table$Left - orf.table$Right ) - 2)/3 )
}



#' @name lorf
#' @title Longest ORF
#' 
#' @description Filtering an ORF-table to keep only the LORFs.
#' 
#' @param orf.table A data.frame listing ORF information, see \code{\link{orfTable}}.
#' 
#' @details For every stop-codon there are usually mutliple possible start-codons in the same reading
#' frame. The LORF (Longest ORF) is defined as the longest of those ORFs sharing the same stop-codon,
#' i.e. the most upstream start-codon matching the stop-codon.
#' 
#' @return An ORF-table with a subset of the rows of the argument \code{orf.table}. In this
#' ORF-table all stop coordinates are unique for each genome sequence.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfTable}}.
#' 
#' @examples Mer her
#' 
#' @export lorf
#' 
lorf <- function( orf.table ){
  ugs <- unique( orf.table$GenomeSequence )
  ot <- orf.table[(orf.table$GenomeSequence == ugs[1]),]
  ot.p <- ot[(ot$Strand > 0),]
  ot.p <- ot.p[order( ot.p$Left ),]
  ot.p <- ot.p[(!duplicated( ot.p$Right )),]
  ot.n <- ot[(ot$Strand < 0),]
  ot.n <- ot.n[order( ot.n$Right, decreasing=T ),]
  ot.n <- ot.n[(!duplicated( ot.n$Left )),]
  os <- rbind( ot.p, ot.n )
  
  if( length( ugs ) > 1 ){
    for( i in 2:length( ugs ) ){
      ot <- orf.table[(orf.table$GenomeSequence == ugs[i]),]
      ot.p <- ot[(ot$Strand > 0),]
      ot.p <- ot.p[order( ot.p$Left ),]
      ot.p <- ot.p[(!duplicated( ot.p$Right )),]
      ot.n <- ot[(ot$Strand < 0),]
      ot.n <- ot.n[order( ot.n$Right, decreasing=T ),]
      ot.n <- ot.n[(!duplicated( ot.n$Left)),]
      os <- rbind( os, ot.p, ot.n )
    }
  }
  return( os )
}




