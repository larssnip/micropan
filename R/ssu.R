#' @name barrnap
#' @title Finding rRNA genes
#' 
#' @description Locating all rRNA genes in genomic DNA using the barrnap software.
#' 
#' @param genome.file A fasta-formatted file with the genome sequence(s).
#' @param bacteria Logical, the genome is either a bacteria (default) or an archea.
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
#' # This example requires the external barrnap software
#' # Using a genome file in this package.
#' xpth <- file.path(path.package("micropan"),"extdata")
#' genome.file <- file.path(xpth,"Example_genome.fasta.xz")
#' 
#' # We need to uncompress it first...
#' tf <- tempfile(fileext=".xz")
#' s <- file.copy(genome.file,tf)
#' tf <- xzuncompress(tf)
#' 
#' # Searching for rRNA sequences, and inspecting
#' gff.table <- barrnap(tf)
#' print(gff.table)
#' 
#' # Retrieving the sequences
#' genome <- readFasta(tf)
#' rRNA.fasta <- gff2fasta(gff.table,genome)
#' 
#' # ...and cleaning...
#' file.remove(tf)
#' }
#' 
#' @export barrnap
#' 
barrnap <- function( genome.file, bacteria=TRUE, cpu=1 ){
  kingdom <- ifelse( bacteria, "bac", "arc" )
  cmd <- paste( "barrnap --quiet --kingdom", kingdom, genome.file, "> barrnap.gff" )
  chr <- NULL
  try( chr <- system( cmd, intern=TRUE), silent=TRUE )
  if( is.null( chr ) ){
    stop( paste('barrnap was not found by R.',
                'Please install barrnap from: https://github.com/tseemann/barrnap',
                'After installation, re-start R and make sure barrnap can be run from R by',
                'the command \'system("barrnap --help")\'.', sep = '\n'))
  } else {
    gff.table <- readGFF( "barrnap.gff" )
    file.remove( "barrnap.gff" )
    return( gff.table )
  }
}
