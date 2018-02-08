#' @name findOrfs
#' @title Finding ORFs in genomes
#' 
#' @description Finds all ORFs in prokaryotic genome sequences.
#' 
#' @param genome A \code{\link{Fasta}} object with the genome sequence(s).
#' @param circular Logical indicating if the genome sequences are completed, circular sequences.
#' 
#' @details A prokaryotic Open Reading Frame (ORF) is defined as a subsequence starting with a  start-codon
#' (ATG, GTG or TTG), followed by an integer number of triplets (codons), and ending with a stop-codon (TAA,
#' TGA or TAG). This function will locate all ORFs in a genome.
#' 
#' The argument \code{genome} will typically have several sequences (chromosomes/plasmids/scaffolds/contigs).
#' It is vital that the \emph{first token} (characters before first space) of every \code{genome$Header} is
#' unique, since this will be used to identify these genome sequences in the output.
#' 
#' Note that for any given stop-codon there are usually multiple start-codons in the same reading
#' frame. This function will return all, i.e. the same stop position may appear multiple times. If
#' you want ORFs with the most upstream start-codon only (LORFs), then filter the output from this function
#' with \code{\link{lorfs}}.
#' 
#' By default the genome sequences are assumed to be linear, i.e. contigs or other incomplete fragments
#' of a genome. In such cases there will usually be some truncated ORFs at each end, i.e. ORFs where either
#' the start- or the stop-codon is lacking. In the \code{gff.table} returned by this function this is marked in the
#' Attributes column. The texts "Truncated=10" or "Truncated=01" indicates truncated at 
#' the Start or End, respectively.
#' 
#' If the supplied \code{genome} is a completed genome, with 
#' circular chromosome/plasmids, set the flag \code{circular=TRUE} and no truncated ORFs will be listed.
#' In cases where an ORF runs across the origin of a circular genome sequences, the Stop coordinate will be
#' larger than the length of the genome sequence. This is in line with the specifications of the GFF3 format, where 
#' a Start cannot be larger than the corresponding End.
#' 
#' @return This function returns a \code{gff.table}, which is simply a \code{data.frame} with columns
#' adhering to the format specified by the GFF3 format, see \code{\link{readGFF}}
#' 
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}, \code{\link{lorfs}}.
#' 
#' @examples
#' \dontrun{
#' # Using a genome file in this package
#' extdata <- file.path(path.package("micropan"),"extdata")
#' genome.file <- "Mpneumoniae_309_genome.fsa"
#' 
#' # We need to uncompress them first...
#' xzuncompress(file.path(extdata,paste(genome.file,".xz",sep="")))
#' 
#' # Reading into R and finding orfs
#' genome <- readFasta(file.path(extdata,genome.file))
#' orf.table <- findOrfs(genome)
#' 
#' # Computing ORF-lengths
#' orf.lengths <- orfLength(orf.table)
#' 
#' # Filtering to retrieve the LORFs only
#' lorf.table <- lorfs(orf.table)
#' lorf.lengths <- orfLength(lorf.table)
#' 
#' # ...and compressing the genome file again...
#' xzcompress(file.path(extdata,genome.file))
#' }
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
  dStrand <- rep( "+", nr )
  dStrand[orf.table$Strand < 0] <- "-"
  dAttribute <- rep( "Truncated=00", nr )
  dAttribute[orf.table$Truncated > 0] <- "Truncated=10"
  dAttribute[orf.table$Truncated < 0] <- "Truncated=01"
  
  gff.table <- data.frame( Seqid=orf.table$Seqid,
                           Source=rep( "micropan::findOrfs", nr ),
                           Type= rep( "ORF", nr ),
                           Start=orf.table$Start,
                           End=orf.table$End,
                           Score=rep( NA, nr ),
                           Strand=dStrand,
                           Phase=rep( 0, nr ),
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
#' @seealso \code{\link{findOrfs}}, \code{\link{lorfs}}.
#' 
#' @examples # See the example in the Help-file for findOrfs.
#' 
#' @export orfLength
#' 
orfLength <- function( gff.table ){
  return( (abs( gff.table$Start - gff.table$End ) - 2)/3 )
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
#' After this filtering the Type variable in \code{gff.table} is changed to \code{"LORF"}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{findOrfs}}.
#' 
#' @examples # See the example in the Help-file for findOrfs.
#' 
#' @export lorfs
#' 
lorfs <- function( gff.table ){
  ugs <- unique( gff.table$Seqid )
  ot <- gff.table[(gff.table$Seqid == ugs[1]),]
  ot.p <- ot[(ot$Strand > 0),]
  ot.p <- ot.p[order( ot.p$Start ),]
  ot.p <- ot.p[(!duplicated( ot.p$End )),]
  ot.n <- ot[(ot$Strand < 0),]
  ot.n <- ot.n[order( ot.n$End, decreasing=T ),]
  ot.n <- ot.n[(!duplicated( ot.n$Start )),]
  os <- rbind( ot.p, ot.n )
  
  if( length( ugs ) > 1 ){
    for( i in 2:length( ugs ) ){
      ot <- gff.table[(gff.table$Seqid == ugs[i]),]
      ot.p <- ot[(ot$Strand > 0),]
      ot.p <- ot.p[order( ot.p$Start ),]
      ot.p <- ot.p[(!duplicated( ot.p$End )),]
      ot.n <- ot[(ot$Strand < 0),]
      ot.n <- ot.n[order( ot.n$End, decreasing=T ),]
      ot.n <- ot.n[(!duplicated( ot.n$Start)),]
      os <- rbind( os, ot.p, ot.n )
    }
  }
  os$Type <- "LORF"
  return( os )
}











# Local function
circularize <- function( ott, NC ){
  tags <- names( NC )
  ugs <- unique( ott$Seqid )
  # otn is the NEW table, where we have spliced the ORFs truncated at each end
  # It is impossible to know how many ORF we will end up with!
  otn <- data.frame( Seqid=NULL, Start=NULL, End=NULL, Strand=NULL, Truncated=NULL, stringsAsFactors=F )
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
        nb1 <- ottg1$End                  # number of bases from origin to start
        ottg2 <- ottg[-ixd,]               # those with truncated stops
        nb2 <- NC[idx] - ottg2$Start + 1   # number of bases before origin
        for( j in 1:length( nb1 ) ){
          idd <- which( ((nb2+nb1[j]) %% 3) == 0 & ottg2$Strand == ottg1$Strand[j] )
          nd <- length( idd )
          if( nd > 0 ){
            otn <- rbind( otn,
                          data.frame( Seqid     = rep( ugs[i], nd ),
                                      Start     = ottg2$Start[idd],
                                      End      = NC[idx]+rep( ottg1$End[j], nd ),
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
