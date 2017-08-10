#' @name orfTable
#' @title Finding ORFs in genomes
#' 
#' @description Finds all ORFs in prokaryotic genome sequences.
#' 
#' @param genome A \code{\link{Fasta}}-object with the genome sequence(s).
#' @param circular Logical indicating if the genome sequences are completed, circular sequences.
#' 
#' @details An Open Reading Frame (ORF) is defined as a subsequence starting with a  start-codon
#' (ATG, GTG or TTG), followed by an integer number of triplets, and ending with a stop-codon (TAA,
#' TGA or TAG). This function will locate all ORFs in a genome, and return a table (\code{data.frame})
#' of information with one row for each ORF, see below. In order to retrieve the sequences, see 
#' \code{\link{orfSequences}}.
#' 
#' Note that for any given stop-codon there are usually multiple start-codons in the same reading
#' frame. This function will return all, i.e. the same stop position may appear multiple times. If
#' you want ORFs with the most upstream start-codon only (LORFs), see \code{\link{lorf}}.
#' 
#' By default the genome sequences are assumed to be linear, i.e. contigs or other incomplete fragments
#' of a genome. In such cases there will usually be som partial ORFs at each end, i.e. ORFs where either
#' the start- or the stop-codon is lacking. If the supplied \code{genome} is a completed genome, with 
#' circular chromosome/plasmids, set the flag \code{circular=TRUE} and no partial ORFs will be listed.
#' You may, however, see some ORFs running across the sequence-start, with the \code{Left}-coordinate
#' larger than the \code{Right}-coordinate.
#' 
#' @return A \code{data.frame} with one row for each ORF and the following columns:
#' \itemize{
#'   \item \code{GenomeSequence}. This is the first token in the \code{genome$Header} texts.
#'   Note that this token must be unique to each \code{genome$Header}-text.
#'   \item \code{Strand}. Either 1 or -1.
#'   \item \code{Left}. The leftmost coordinate. This is the start if the ORF is on the 1 strand, but
#'   the stop if it is on the -1 strand.
#'   \item \code{Right}. The rightmost coordinate. This is the stop if the ORF is on the 1 strand, but
#'   the start if it is on the -1 strand.
#'   \item \code{Partial}. This is 1 if the ORF in incomplete at the start of the genome sequence and -1
#'   it is incomplete at the end. All complete ORFs have value 0.
#' }
#' 
#' Note that the \code{Left} and \code{Right} coordinates specify the leftmost and rightmost base in the
#' ORF regardless of strand. If \code{Strand=1} the start-codon is in the positions (\code{Left}, 
#' \code{Left}+1,\code{Left}+2) and the stop-codon in (\code{Right}-2, \code{Right}-1, \code{Right}).
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfSequences}}, \code{\link{lorf}}.
#' 
#' @examples Mer her
#' 
#' @useDynLib micropan
#' @importFrom Rcpp evalCpp
#' 
#' @export orfTable
#' 
orfTable <- function( genome, circular=F ){
  tags <- sapply( strsplit( genome$Header, split=" " ), function(x){x[1]} )
  if( length( unique( tags ) ) != length( tags ) ) stop( "First token in the Headers must be unique!" )
  orf.table <- ORF_index( tags, genome$Sequence )

  if( circular ){
    is.part <- (orf.table$Partial != 0)
    otp <- orf.table[is.part,]
    ugs <- as.character( unique( otp$GenomeSequence ) )
    otn <- data.frame( GenomeSequence=NULL, Strand=NULL, Left=NULL, Right=NULL, Partial=NULL, stringsAsFactors=F )
    for( i in 1:length(ugs) ){
      NC <- nchar( genome$Sequence[grep( ugs[i], genome$Header, fixed=T )] )
      otpg <- otp[which( otp$GenomeSequence == ugs[i] ),]
      nr <- nrow( otpg )
      if( nr > 0 ){
        ixd <- which( otpg$Partial > 0 )   # incomplete start
        ni <- length( ixd )
        if( (ni > 0) & (ni < nr) ){
          otpg1 <- otpg[ixd,]
          nb1 <- otpg1$Right         # number of bases after start
          otpg2 <- otpg[-ixd,]
          nb2 <- NC - otpg2$Left + 1 # number of bases before start
          for( j in 1:length( nb1 ) ){
            idd <- which( ((nb2+nb1[j]) %% 3) == 0 & otpg2$Strand == otpg1$Strand[j] )
            nd <- length( idd )
            if( nd > 0 ){
              otn <- rbind( otn,
                            data.frame( GenomeSequence= rep( ugs[i], nd ),
                                        Strand        = otpg2$Strand[idd],
                                        Left          = otpg2$Left[idd],
                                        Right         = rep( otpg1$Right[j], nd ),
                                        Partial       = rep( 0, nd ),
                                        stringsAsFactors=F ) )
            }
          }
        }
      }
    }
    orf.table <- rbind( orf.table[!is.part,], otn )
  }
  return( orf.table )
}



#' @name orfSequences
#' @title Retrieving ORFs
#' 
#' @description Retrieving ORF sequences from a genome.
#' 
#' @param orf.table A data.frame with ORF information, see \code{\link{orfTable}}.
#' @param genome A \code{\link{Fasta}}-object with the genome sequence(s).
#' 
#' @details For each ORF listed in \code{orf.table} the sequence is retrieved from the
#' genome. Every \code{GenomeSequence} tag in the \code{orf.table} must match the first token
#' in one of the \code{genome$Header} texts.
#' 
#' @return A \code{\link{Fasta}} object with one row for each ORF. The \code{Header} for each
#' ORF is a summary of the \code{orf.table} information for the same ORF.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{orfTable}}, \code{\link{lorf}}.
#' 
#' @examples Mer her
#' 
#' @useDynLib micropan
#' @importFrom Rcpp evalCpp
#' 
#' @export orfSequences
#' 
orfSequences <- function( orf.table, genome ){
  utag <- unique( orf.table$GenomeSequence )
  tagz <- unique( sapply( strsplit( genome$Header, split=" " ), function(x){x[1]} ) )
  if( sum( is.na( match( utag, tagz ) ) ) > 0 )
    stop( "GenomeSequence tags in the orf.table do not match the genome$Header" )
  
  n <- nrow( orf.table )
  fobj <- data.frame( Header=orfSignature( orf.table ), 
                      Sequence=rep( "", n ),
                      stringsAsFactors=F )
  for( i in 1:length( tagz ) ){
    idx <- which( orf.table$GenomeSequence == tagz[i] )
    seq <- extractSeq( genome$Sequence[i], 
                      orf.table$Left[idx], 
                      orf.table$Right[idx], 
                      orf.table$Strand[idx] )
    idd <- which( orf.table$Strand[idx] < 0 )
    seq[idd] <- reverseComplement( seq[idd] )
    fobj$Sequence[idx] <- seq
  }
  class( fobj ) <- c( "Fasta", "data.frame" )
  return( fobj )
}

#' @name orfSignature
#' @title ORF signature text
#' 
#' @description Making a signature text from ORF data.
#' 
#' @param orf.table A data.frame with ORF information, see \code{\link{orfTable}}.
#' 
#' @details For each ORF listed in \code{orf.table} a text is created by pasting these
#' data together, adding some explanatory text. This function is used by \code{link{orfSequences}}
#' to create the header-lines for the \code{\link{Fasta}} object when retrieving the ORF sequences.
#' 
#' @return A vector of texts, one for each row in \code{orf.table}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{orfTable}}, \code{\link{orfSequences}}.
#' 
#' @examples Mer her
#' 
#' @useDynLib micropan
#' 
#' @export orfSignature
#' 
orfSignature <- function( orf.table ){
  desc <- paste( "GenomeSequence=", orf.table$GenomeSequence,
                 ";Strand=", orf.table$Strand,
                 ";Left=", orf.table$Left,
                 ";Right=", orf.table$Right, sep="" )
  if( "Partial" %in% colnames( orf.table ) ){
    desc <- paste( desc, ";Partial=", orf.table$Partial, sep="" )
  }
  return( desc )
}
