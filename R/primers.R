#' @name amplicon
#' @title Primer matching
#' 
#' @description Extracts subsequences from DNA that matches a specified primer pair.
#' 
#' @param dna Character vector containing the DNA sequences.
#' @param forward String specifying the forward primer.
#' @param reverse String specifying the reverse primer.
#' 
#' @details An amplicon is a subsequence limited by a matching pair of short oligos, called primers.
#' 
#' The \code{forward} primer is a short DNA sequence in the 5' to 3' direction. This can match on both strands of
#' the \code{dna} sequence. The \code{reverse} primer is also a short DNA sequence in 5' to 3' direction, and
#' can also match on both strands of \code{dna}.
#' 
#' For a \code{dna} sequence to produce an amplicon there must be an exact match of the \code{forward}
#' on one strand followed by an exact match of the \code{reverse} on the other strand. The amplicon is
#' the subsequence starting with the \code{forward} and ending with the \code{reverse} primer.
#' 
#' Both primers may contain ambiguity symbols according to the IUPAC standard. 
#' Primers are matched by \code{\link{gregexpr}}, which will not register self-overlapping matches. 
#' In case of multiple (non-overlapping) matches, this function will return all possible amplicons resulting 
#' from the primer matching.
#' 
#' @return A list with the same number of elements as the argument \code{dna}. Each list element
#' contains a string vector with all amplicons resulting from the primer matching. If there is no
#' primer pair match the corresponding string vector is empty.
#' 
#' @author Lars Snipen.
#' 
#' @examples
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fasta")
#' fdta <- readFasta(ex.file)
#' amp.lst.1 <- amplicon( fdta$Sequence, forward="AAATTC", reverse="CCAGA" )
#' amp.lst.2 <- amplicon( fdta$Sequence, forward="AANNTC", reverse="CCNGT" ) # more matches due to N's
#' 
#' @export amplicon
#' 
amplicon <- function( dna, forward, reverse ){
  NC <- nchar( dna )
  rvdna <- reverseComplement( dna )
  fwd <- iupac2regex( forward )
  rvs <- iupac2regex( reverse )
  
  ff <- gregexpr( fwd, dna )
  rr <- gregexpr( rvs, rvdna )
  fr <- gregexpr( fwd, rvdna )
  rf <- gregexpr( rvs, dna )
  amplicons <- lapply( 1:length(dna), function(i){
    amps <- character( 0 )
    ss <- ff[[i]]
    ee <- rr[[i]]
    if( ss[1]>0 & ee[1]>0 ){
      ee <- NC[i] - ee + 1
      se <- expand.grid( ss, ee )
      se <- se[which( se[,1]<se[,2] ),]
      if( nrow( se ) >0 ){
        amps <- c( amps,
                   substring( dna[i], se[,1], se[,2] ) )
      }
    }
    ss <- fr[[i]]
    ee <- rf[[i]]
    if( ss[1]>0 & ee[1]>0 ){
      ee <- NC[i] - ee + 1
      se <- expand.grid( ss, ee )
      se <- se[which( se[,1]<se[,2] ),]
      if( nrow( se ) > 0 ){
        amps <- c( amps,
                   reverseComplement( substring( rvdna[i], se[,1], se[,2] ) ) )
      }
    }
    return( amps )
  })
  return( amplicons )
}



