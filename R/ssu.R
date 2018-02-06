#' @name barrnap
#' @title Finding rRNA genes
#' 
#' @description Locating all rRNA genes in genomic DNA using the barrnap software.
#' 
#' @param genome.file A fasta-formatted file with the genome sequence(s).
#' @param bacteria. Logical, the genome is either a bacteria (default) or an archea.
#' @param cpu Number of CPUs to use, default is 1.
#' 
#' @details The external software barrnap is used to scan through a prokaryotic genome to detect the
#' rRNA genes (5S, 16S, 23S). This free software can be installed from https://github.com/tseemann/barrnap.
#' 
#' @return A \code{gff.table} (see \code{\link{readGFF}} for details) with one row for each detected
#' rRNA sequence. 
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}.
#' 
#' @examples
#' # Using a genome FASTA file in this package.
#' # We need to uncompress it first...
#' \dontrun{
#' extdata.path <- file.path(path.package("micropan"),"extdata")
#' genome.file <- xzuncompress(file.path(extdata.path, "Mpneumoniae_309_genome.fsa.xz"))
#' 
#' #...searching for rRNA sequences...
#' ssu.table <- barrnap( genome.file )
#' 
#' # Inspecting the hits
#' print(ssu.table)
#' 
#' # Retrieving the sequences
#' genome <- readFasta(genome.file)
#' rRNA.fasta <- gff2fasta(ssu.table, genome)
#' 
#' # ...and compressing the genome FASTA file again...
#' genome.file.xz <- xzcompress(genome.file)
#' }
#' 
#' @export barrnap
#' 
barrnap <- function( genome.file, bacteria=TRUE, cpu=1 ){
  if( available.external( "barrnap" ) ){
    kingdom <- ifelse( Bacteria, "bac", "arc" )
    cmd <- paste( "barrnap --quiet --kingdom", kingdom, genome.file, "> barrnap.gff" )
    system( cmd )
    gff.table <- readGFF( "barrnap.gff" )
    return( gff.table )
  }
}
