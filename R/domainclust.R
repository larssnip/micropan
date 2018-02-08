#' @name dClust
#' @title Clustering sequences based on domain sequence
#' 
#' @description Proteins are clustered by their sequence of protein domains. A domain sequence is the
#' ordered sequence of domains in the protein. All proteins having identical domain sequence are assigned
#' to the same cluster.
#' 
#' @param hmmer.table A \code{data.frame} of results from a \code{\link{hmmerScan}} against a domain database.
#' 
#' @details A domain sequence is simply the ordered list of domains occurring in a protein. Not all proteins
#' contain known domains, but those who do will have from one to several domains, and these can be ordered
#' forming a sequence. Since domains can be more or less conserved, two proteins can be quite different in
#' their amino acid sequence, and still share the same domains. Describing, and grouping, proteins by their
#' domain sequence was proposed by Snipen & Ussery (2012) as an alternative to clusters based on pairwise
#' alignments, see \code{\link{bClust}}. Domain sequence clusters are less influenced by gene prediction errors.
#' 
#' The input is a \code{data.frame} of the type produced by \code{\link{readHmmer}}. Typically, it is the
#' result of scanning proteins (using \code{\link{hmmerScan}}) against Pfam-A or any other HMMER3 database
#' of protein domains. It is highly reccomended that you remove overlapping hits in \samp{hmmer.table} before
#' you pass it as input to \code{\link{dClust}}. Use the function \code{\link{hmmerCleanOverlap}} for this.
#' Overlapping hits are in some cases real hits, but often the poorest of them are artifacts.
#' 
#' @return The output is a numeric vector with one element for each unique sequence in the \samp{Query}
#' column of the input \samp{hmmer.table}. Sequences with identical number belong to the same cluster. The
#' name of each element identifies the sequence.
#' 
#' This vector also has an attribute called \samp{cluster.info} which is a character vector containing the
#' domain sequences. The first element is the domain sequence for cluster 1, the second for cluster 2, etc.
#' In this way you can, in addition to clustering the sequences, also see which domains the sequences of a
#' particular cluster share.
#' 
#' @references Snipen, L. Ussery, D.W. (2012). A domain sequence approach to pangenomics: Applications
#' to Escherichia coli. F1000 Research, 1:19.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panPrep}}, \code{\link{hmmerScan}}, \code{\link{readHmmer}},
#' \code{\link{hmmerCleanOverlap}}, \code{\link{bClust}}.
#' 
#' @examples 
#' \dontrun{
#' # Using HMMER3 result files in the micropan package
#' extdata <- file.path(path.package("micropan"),"extdata")
#' hmm.files <- c("GID1_vs_Pfam-A.hmm.txt","GID2_vs_Pfam-A.hmm.txt","GID3_vs_Pfam-A.hmm.txt")
#' 
#' # We need to uncompress them first...
#' pth <- lapply(file.path(extdata,paste(hmm.files,".xz",sep="")),xzuncompress)
#' 
#' # Reading the HMMER3 results, cleaning overlaps...
#' hmmer.table <- NULL
#' for(i in 1:3){
#'   htab <- readHmmer(file.path(extdata,hmm.files[i]))
#'   htab <- hmmerCleanOverlap(htab)
#'   hmmer.table <- rbind(hmmer.table,htab)
#' }
#' 
#' # The clustering
#' clustering.domains <- dClust(hmmer.table)
#' 
#' # ...and compressing the result files again...
#' pth <- lapply(file.path(extdata,hmm.files),xzcompress)
#' }
#' 
#' @export dClust
#' 
dClust <- function( hmmer.table ){
  cat( "dClust:\n" )
  cat( "...hmmer.table contains", length( unique( hmmer.table$Query ) ), "proteins...\n" )
  cat( "...with hits against", length( unique( hmmer.table$Hit ) ), "HMMs...\n" )
  hmmer.table <- hmmer.table[order( hmmer.table$Start ),]
  dseq <- sort( unlist( tapply( hmmer.table$Hit, hmmer.table$Query, function( x ){ paste( x, collapse="," ) } ) ) )
  seq <- names( dseq )
  dsc <- as.numeric( factor( as.vector( dseq ) ) )
  names( dsc ) <- seq
  cat( "...ended with", length( unique( dsc ) ), "clusters, largest cluster has", max( table( dsc ) ), "members\n" )
  attr( dsc, "cluster.info" ) <- unique( as.vector( dseq ) )
  return( dsc )
}


#' @name hmmerCleanOverlap
#' @title Removing overlapping hits from HMMER3 scans
#' 
#' @description Removing hits to avoid overlapping HMMs on the same protein sequence.
#' 
#' @param hmmer.table A \code{data.frame} with \code{\link{hmmerScan}} results, see \code{\link{readHmmer}}.
#' 
#' @details  When scanning sequences against a profile HMM database using \code{\link{hmmerScan}}, we
#' often find that several patterns (HMMs) match in the same region of the query sequence, i.e. we have
#' overlapping hits. The function \code{\link{hmmerCleanOverlap}} will remove the poorest overlapping hit
#' in a recursive way such that all overlaps are eliminated.
#' 
#' The input is a \code{data.frame} of the type produced by \code{\link{readHmmer}}.
#' 
#' @return A \code{data.frame} which is a subset of the input, where some rows have been deleted to
#' avoid overlapping hits.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{hmmerScan}}, \code{\link{readHmmer}}, \code{\link{dClust}}.
#' 
#' @examples # See the example in the Help-file for dClust.
#' 
#' @export
hmmerCleanOverlap <- function( hmmer.table ){
  qt <- table( hmmer.table$Query )
  cat( "There are", length( qt ), "proteins in this hmmer.table...\n" )
  hmmer.table <- hmmer.table[order( hmmer.table$Start ),]
  if( max(qt) > 1 ){
    idx <- which( qt > 1 )
    multi <- names( qt )[idx]
    nm <- length( multi )
    cat( "There are", nm, "proteins with multiple hits, resolving overlaps:\n")
    keep <- rep( T, dim( hmmer.table )[1] )
    for( i in 1:nm ){
      idx <- which( hmmer.table$Query == multi[i] )
      keep[idx] <- nonoverlap( hmmer.table[idx,] )
      if( (i/100) == round(i/100) ) cat( "." )
    }
    hmmer.table <- hmmer.table[keep,]
  }
  return( hmmer.table )
}

nonoverlap <- function( hmmer.table ){
  nh <- dim( hmmer.table )[1]
  keep <- rep( T, nh )
  ht <- hmmer.table[keep,]
  dif <- ht$Start[2:nh] - ht$Stop[1:(nh-1)]
  ido <- which( dif <= 0 )
  while( (nh > 1) & (length( ido ) > 0) ){
    idx <- unique( c( ido, ido+1 ) )
    idd <- which( ht$Evalue[idx] == max( ht$Evalue[idx] ) )
    map <- which( keep )
    keep[map[idx[idd[1]]]] <- F
    ht <- hmmer.table[keep,]
    nh <- dim( ht )[1]
    if( nh > 1 ){
      dif <- ht$Start[2:nh] - ht$Stop[1:(nh-1)]
      ido <- which( dif <= 0 )
    }
  }
  return( keep )
}












