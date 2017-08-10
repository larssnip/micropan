#' @name cmsearch16S
#' @title Finding 16S sequences
#' 
#' @description Locating all occurrences of the 16S sequence in genomic DNA
#' using the Infernal software.
#' 
#' @param genome.file A fasta-formatted file with the genome sequence(s).
#' @param cpu Number of CPUs to use, default is 1.
#' @param spacer Minimum number of bases between two distinct 16S genes on the same genome sequence.
#' @param verbose Logical, if \code{TRUE} progress is reported during computations (default=TRUE).
#' 
#' @details The software Infernal (Nawrocki & Eddy, 2013) must be installed and available on the system. Test 
#' this by typing \code{system("cmsearch -h")} in the Console, and some sensible output should be produced. 
#' For more details on Infernal, see http://eddylab.org/infernal/.
#' 
#' The module \code{cmsearch} of the Infernal software is used to search with the correlation
#' models for both bacteria and archaea from the Rfam database (http://rfam.xfam.org/) against the specified
#' genome. All significant hits (E-value below 1.0) are considered, but only unique hits are reported. In case
#' of multiple hits in the same genome sequence, hits are clustered by their location, using hirerachical 
#' clustering with single linkage. Thus, if hits are separated by less than \code{spacer} nucleotides, they are
#' clustered, and only the hit with the largest score in each cluster is reported. Setting \code{spacer=0} will
#' result in all hits being reported.
#' 
#' A full-length 16S sequence should be at least 
#' 1200 bases, but (much) shorter hits may be reported. Incomplete genomes or problematic
#' assemblies will typically produce these shorter hits.
#' 
#' @return A \code{data.frame} with one row for each detected 16S sequence. Similar to ORFs
#' this table specifies the \code{GenomeSequence}, \code{Strand},
#' \code{Left} and \code{Right} for each sequence, see \code{\link{orfTable}}
#' for details on this. There is no \code{Partial} information for 16S sequences, as their
#' exact start and stop varies too much.
#' 
#' The output from this function can be used, together with a \code{\link{Fasta}} object with the genome
#' to retrieve the 16S sequences by \code{\link{orfSequences}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @references E.P. Nawrocki and S.R. Eddy,  Infernal 1.1: 100-fold faster RNA homology searches, 
#' Bioinformatics 29:2933-2935 (2013). 
#' 
#' @seealso \code{\link{orfTable}}, \code{\link{orfSequences}}.
#' 
#' @examples
#' # Using a FASTA file in this package.
#' # We need to uncompress it first...
#' \dontrun{
#' extdata.path <- file.path(path.package("micropan"),"extdata")
#' genome.file <- xzuncompress(file.path(extdata.path, "Mpneumoniae_309_genome.fsa.xz"))
#' #...searching for 16S sequences...
#' ssu.table <- cmsearch16S( genome.file )
#' 
#' # Inspecting the hits
#' print(ssu.table)
#' 
#' # Retrieving the sequences
#' gdta <- readFasta(genome.file)
#' fdta.16S <- orfSequences(ssu.table, gdta)
#' 
#' # ...and compressing the FASTA file again...
#' genome.file.xz <- xzcompress(genome.file)
#' }
#' 
#' @export cmsearch16S
#' 
cmsearch16S <- function( genome.file, cpu=1, spacer=2000, verbose=TRUE ){
  if( available.external( "infernal" ) ){
    if( verbose ) cat( "cmsearch16S: " )
    cm.file <- file.path( path.package( "micropan" ), "extdata/ssu.cm" )
    cmd <- paste( "cmsearch -E 1.0 --tblout res.txt --noali -o logfile.txt --cpu", cpu, 
                  cm.file, genome.file )
    system( cmd )
    lines <- readLines( "res.txt" )
    file.remove( "res.txt", "logfile.txt" )
    idx <- grep( "^#", lines, invert=T )
    if( length( idx ) == 0 ){
      if( verbose ) cat( "  found no hits!" )
      ssu.tab <- data.frame( Genome.sequence=NULL,
                             Left=NULL, Right=NULL, Strand=NULL,
                             stringsAsFactors=F )
    } else {
      rlines <- lines[idx]
      lst <- strsplit( rlines, split=" +" )
      h <- sapply( lst, function(x){x[1]} )
      q <- sapply( lst, function(x){x[3]} )
      hs <- as.integer( sapply( lst, function(x){x[8]} ) )
      he <- as.integer( sapply( lst, function(x){x[9]} ) )
      st <- sapply( lst, function(x){x[10]} )
      sc <- as.numeric( sapply( lst, function(x){x[15]} ) )
      rtab <- data.frame( Hit=h, Query=q, Hit.start=hs, Hit.end=he,
                          Strand=st, Score=sc,
                          stringsAsFactors=F )
      score.mean <- sort( tapply( rtab$Score, rtab$Query, mean ), decreasing=T )
      dom <- names( score.mean )[1]
      rtab <- rtab[rtab$Query==dom,]
      if( verbose ){
        cat( "  found", nrow( rtab ), "hits..." )
        if( rtab$Query[1] == "SSU_rRNA_bacteria" ) cat( "(looks like a bacteria)" )
        else cat( "(looks like an archaea)" )
      }
      uhits <- unique( rtab$Hit )
      ssu.tab <- data.frame( GenomeSequence=NULL, Strand=NULL, Left=NULL, Right=NULL, stringsAsFactors=F )
      for( i in 1:length( uhits ) ){
        idx <- which( rtab$Hit == uhits[i] )
        rt <- rtab[idx,]
        if( nrow( rt ) > 1 ){ # multiple hits, needs clustering to resolve unique hits
          d1 <- t( sapply( rt$Hit.start,function(x){pmin(abs(x-rt$Hit.start),abs(x-rt$Hit.end))}) )
          d2 <- t( sapply( rt$Hit.end,  function(x){pmin(abs(x-rt$Hit.start),abs(x-rt$Hit.end))}) )
          d <- pmin( d1, d2 )
          tree <- hclust( as.dist( d ), method="single" )
          clst <- cutree( tree, h=spacer )
          uclst <- unique( clst )
          nnu <- length( uclst )
          ssu <- data.frame( GenomeSequence=rep("",nnu),
                             Strand=rep(0,nnu),
                             Left=rep(0,nnu),
                             Right=rep(0,nnu),
                             stringsAsFactors=F )
          for( j in 1:nnu ){
            idx <- which( clst == uclst[j] )
            idd <- which( rt$Score[idx] == max( rt$Score[idx] ) )
            ssu$GenomeSequence[j] <- rt$Hit[idx[idd[1]]]
            ssu$Strand[j] <- ifelse( rt$Strand[idx[idd[1]]] == "+", 1, -1 )
            ssu$Left[j]  <- ifelse( ssu$Strand[j] == 1, rt$Hit.start[idx[idd[1]]], rt$Hit.end[idx[idd[1]]] )
            ssu$Right[j] <- ifelse( ssu$Strand[j] == 1, rt$Hit.end[idx[idd[1]]], rt$Hit.start[idx[idd[1]]] )
          }
        } else { # single hit only
          ssu <- data.frame( GenomeSequence=rt$Hit,
                             Strand=ifelse( rt$Strand == "+", 1, -1 ),
                             Left=ifelse( rt$Strand == "+", rt$Hit.start, rt$Hit.end ),
                             Right=ifelse( rt$Strand == "+", rt$Hit.end, rt$Hit.start ), 
                             stringsAsFactors=F )
        }
        ssu.tab <- rbind( ssu.tab, ssu )
      } # new uhit
    }
    if( verbose ) cat( "...ended with", nrow( ssu.tab ), "unique hits\n" )
    return( ssu.tab )
  }
}
  
  
  
## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function( what ){
  if(what == "infernal"){
    chr <- NULL
    try(chr <- system('cmsearch -h', intern = TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('Infernal was not found by R.',
                 'Please install Infernal from: http://eddylab.org/infernal/.',
                 'After installation, make sure Infernal can be run from a shell/terminal, ',
                 'using the command \'cmsearch\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}
  