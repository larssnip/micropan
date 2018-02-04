#' @name findOrfs
#' @title Finding ORFs in genomes
#' 
#' @description Finds all ORFs in prokaryotic genome sequences.
#' 
#' @param genome A \code{\link{Fasta}} object with the genome sequence(s).
#' @param circular Logical indicating if the genome sequences are completed, circular sequences.
#' 
#' @details An Open Reading Frame (ORF) is defined as a subsequence starting with a  start-codon
#' (ATG, GTG or TTG), followed by an integer number of triplets, and ending with a stop-codon (TAA,
#' TGA or TAG). This function will locate all ORFs in a genome, and returns a \code{gff.table}, which is
#' explained below.
#' 
#' The argument \code{genome} will typically have several sequences (chromosomes/plasmids/scaffolds/contigs).
#' It is vital that the \emph{first token} (characters before first space) of every \code{genome$Header} is
#' unique, since this will be used to identify these genome sequences in the output.
#' 
#' Note that for any given stop-codon there are usually multiple start-codons in the same reading
#' frame. This function will return all, i.e. the same stop position may appear multiple times. If
#' you want ORFs with the most upstream start-codon only (LORFs), see \code{\link{lorfs}}.
#' 
#' By default the genome sequences are assumed to be linear, i.e. contigs or other incomplete fragments
#' of a genome. In such cases there will usually be some truncated ORFs at each end, i.e. ORFs where either
#' the start- or the stop-codon is lacking. If the supplied \code{genome} is a completed genome, with 
#' circular chromosome/plasmids, set the flag \code{circular=TRUE} and no truncated ORFs will be listed.
#' In cases where an ORF runs across the origin of a circular genome sequences, the Stop coordinate will be
#' larger than the length of the genome sequence. This is in line with the specifications of the GFF3 format, where 
#' a Start cannot be larger than the corresponding Stop
#' 
#' @return This function returns a \code{gff.table}, which is simply a \code{data.frame} with columns
#' adhering, more or less, to the format specified by the GFF3 format, see
#' https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details. There is
#' one row for each ORF.
#' 
#' The following columns are in the \code{gff.table} returned by \code{findOrfs}:
#' \itemize{
#'   \item \code{Seqid}. This is the first token in the \code{genome$Header} texts, i.e. an identifier for
#'   the genome sequence/scaffold/contig where the ORF is found. Note that this implies the first token
#'   in all \code{genome$Header} texts must be unique!
#'   \item \code{Type} All features are of type \code{"ORF"}.
#'   \item \code{Start}. The leftmost coordinate. This is the start if the ORF is on the Sense strand, but
#'   the stop if it is on the Antisense strand.
#'   \item \code{Stop}. The rightmost coordinate. This is the stop if the ORF is on the Sense strand, but
#'   the start if it is on the Antisense strand.
#'   \item \code{Strand}. A \code{"+"} indicates Sense strand, a \code{"-"} Antisense.
#'   \item \code{Attributes}. A text indicating if the ORF is truncated. If \code{"truncated=00"} the ORF
#'   is complete. If \code{"truncated=10"} it is 
#' }
#' These are columns 1,3,4,5,7 and 9 of the GFF3 format. Columns 2 (source), 6 (score) and 8 (phase) are not 
#' included, since they will not contain informative data from this function. 
#' 
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{gff2fasta}}, \code{\link{lorfs}}.
#' 
#' @examples
#' # Reading a small genome into R (uncompressing it first):
#' extdata.path <- file.path(path.package("micropan"),"extdata")
#' pth <- lapply( file.path( extdata.path, "Mpneumoniae_309_genome.fsa.xz" ), xzuncompress )  # uncompressing file...
#' genome <- readFasta( file.path( extdata.path, "Mpneumoniae_309_genome.fsa" ) )
#' pth <- lapply( file.path( extdata.path, "Mpneumoniae_309_genome.fsa" ), xzcompress )   # ...and compressing it again...
#' 
#' # Finding ORFs
#' gff.tab <- findOrfs( genome )
#' 
#' @useDynLib micropan
#' @importFrom Rcpp evalCpp
#' 
#' @export findOrfs
#' 
findOrfs <- function( genome, circular=F ){
  tags <- sapply( strsplit( genome$Header, split=" " ), function(x){x[1]} )
  if( length( unique( tags ) ) != length( tags ) ) stop( "First token in the Headers must be unique!" )
  NC <- nchar( genome$Sequence )
  names( NC ) <- tags
  orf.table <- ORF_index( tags, genome$Sequence )
  orf.table$Seqid <- as.character( orf.table$Seqid )
  
  if( circular ){
    idx <- which( orf.table$Truncated != 0 )
    otn <- circularize( orf.table[idx,], NC )
    orf.table <- rbind( orf.table[-idx,],
                        otn )
  }
  
  nr <- nrow( orf.table )
  dType <- rep( "ORF", nr )
  dStrand <- rep( "+", nr )
  dStrand[orf.table$Strand < 0] <- "-"
  dAttribute <- rep( "Truncated=00", nr )
  dAttribute[orf.table$Truncated > 0] <- "Truncated=10"
  dAttribute[orf.table$Truncated < 0] <- "Truncated=01"
  
  gff.table <- data.frame( Seqid=orf.table$Seqid,
                           Type=dType,
                           Start=orf.table$Start,
                           Stop=orf.table$Stop,
                           Strand=dStrand,
                           Attributes=dAttribute,
                           stringsAsFactors=F )
  return( gff.table )
}



#' @name orfLength
#' @title ORF lengths
#' 
#' @description Computes the lengths of all ORFs in a \code{gff.table}.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' 
#' @details See \code{\link{findOrfs}} for more on \code{gff.table}s.
#' 
#' @return A vector of ORF lengths, measured as the number of amino acids after translation.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{findOrfs}}.
#' 
#' @examples Mer her
#' 
#' @export orfLength
#' 
orfLength <- function( gff.table ){
  return( (abs( gff.table$Start - gff.table$Stop ) - 2)/3 )
}



#' @name lorfs
#' @title Longest ORF
#' 
#' @description Filtering a \code{gff.table} with ORF information to keep only the LORFs.
#' 
#' @param gff.table A \code{gff.table} (\code{data.frame}) with genomic features information.
#' 
#' @details For every stop-codon there are usually mutliple possible start-codons in the same reading
#' frame (nested ORFs). The LORF (Longest ORF) is defined as the longest of these nested ORFs,
#' i.e. the ORF starting at the most upstream start-codon matching the stop-codon.
#' 
#' @return A \code{gff.table} with a subset of the rows of the argument \code{gff.table}. 
#' After this filtering all stop coordinates are unique for each genome sequence. Also, the 
#' \code{Type} variable is changed to \code{"LORF"}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{findOrfs}}.
#' 
#' @examples Mer her
#' 
#' @export lorfs
#' 
lorfs <- function( gff.table ){
  ugs <- unique( gff.table$Seqid )
  ot <- gff.table[(gff.table$Seqid == ugs[1]),]
  ot.p <- ot[(ot$Strand > 0),]
  ot.p <- ot.p[order( ot.p$Start ),]
  ot.p <- ot.p[(!duplicated( ot.p$Stop )),]
  ot.n <- ot[(ot$Strand < 0),]
  ot.n <- ot.n[order( ot.n$Stop, decreasing=T ),]
  ot.n <- ot.n[(!duplicated( ot.n$Start )),]
  os <- rbind( ot.p, ot.n )
  
  if( length( ugs ) > 1 ){
    for( i in 2:length( ugs ) ){
      ot <- gff.table[(gff.table$Seqid == ugs[i]),]
      ot.p <- ot[(ot$Strand > 0),]
      ot.p <- ot.p[order( ot.p$Start ),]
      ot.p <- ot.p[(!duplicated( ot.p$Stop )),]
      ot.n <- ot[(ot$Strand < 0),]
      ot.n <- ot.n[order( ot.n$Stop, decreasing=T ),]
      ot.n <- ot.n[(!duplicated( ot.n$Start)),]
      os <- rbind( os, ot.p, ot.n )
    }
  }
  os$Type <- "LORF"
  return( os )
}












circularize <- function( ott, NC ){
  tags <- names( NC )
  ugs <- unique( ott$Seqid )
  # otn is the NEW table, where we have spliced the ORFs truncated at each end
  # It is impossible to know how many ORF we will end up with!
  otn <- data.frame( Seqid=NULL, Start=NULL, Stop=NULL, Strand=NULL, Truncated=NULL, stringsAsFactors=F )
  for( i in 1:length(ugs) ){
    idx <- which( tags == ugs[i] )                # ugs[i] is genome sequence idx
    ottg <- ott[which( ott$Seqid == ugs[i] ),]    # orfs from sequence idx
    nr <- nrow( ottg )
    if( nr > 0 ){                          # there are truncated ORFs for this genome sequence
      ixd <- which( ottg$Truncated > 0 )   # incomplete starts
      ni <- length( ixd )
      if( (ni > 0) & (ni < nr) ){          # we need some with truncated starts, but also some
                                           # with truncated stops!
        ottg1 <- ottg[ixd,]                # those with truncated starts
        nb1 <- ottg1$Stop                  # number of bases from origin to start
        ottg2 <- ottg[-ixd,]               # those with truncated stops
        nb2 <- NC[idx] - ottg2$Start + 1   # number of bases before origin
        for( j in 1:length( nb1 ) ){
          idd <- which( ((nb2+nb1[j]) %% 3) == 0 & ottg2$Strand == ottg1$Strand[j] )
          nd <- length( idd )
          if( nd > 0 ){
            otn <- rbind( otn,
                          data.frame( Seqid     = rep( ugs[i], nd ),
                                      Start     = ottg2$Start[idd],
                                      Stop      = NC[idx]+rep( ottg1$Stop[j], nd ),
                                      Strand    = ottg2$Strand[idd],
                                      Truncated = rep( 0, nd ),
                                      stringsAsFactors=F ) )
          }
        }
      }
    }
  }
  
  
  return( otn )
}
