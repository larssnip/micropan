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
#' @note The barrnap software must be installed on the system for this function to work, i.e. the command
#' \samp{system("barrnap --help")} must be recognized as a valid command if you run it in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}, \code{\link{gff2fasta}}.
#' 
#' @examples
#' \dontrun{
#' # Using a genome file in this package.
#' extdata <- file.path(path.package("micropan"),"extdata")
#' genome.file <- "Mpneumoniae_309_genome.fsa"
#' 
#' # We need to uncompress it first...
#' xzuncompress(file.path(extdata,paste(genome.file,".xz",sep="")))
#' 
#' # Searching for rRNA sequences, and inspecting
#' gff.table <- barrnap(file.path(extdata,genome.file))
#' print(gff.table)
#' 
#' # Retrieving the sequences
#' genome <- readFasta(file.path(extdata,genome.file))
#' rRNA.fasta <- gff2fasta(ssu.table,genome)
#' 
#' # ...and compressing the genome FASTA file again...
#' xzcompress(file.path(extdata,genome.file))
#' }
#' 
#' @export barrnap
#' 
barrnap <- function( genome.file, bacteria=TRUE, cpu=1 ){
  if( available.external( "barrnap" ) ){
    kingdom <- ifelse( bacteria, "bac", "arc" )
    cmd <- paste( "barrnap --quiet --kingdom", kingdom, genome.file, "> barrnap.gff" )
    system( cmd )
    gff.table <- readGFF( "barrnap.gff" )
    return( gff.table )
  }
}
