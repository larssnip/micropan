#' @name prodigalPredict
#' @title Gene predictions using Prodigal
#' 
#' @description Finds coding genes in a genome, using the Prodigal software, and outputs them as a FASTA file.
#' 
#' @param genome.file Name of a FASTA formatted file with all the DNA sequences for a genome (chromosomes,
#' plasmids, contigs etc.).
#' @param prot.file Name of output file. Predicted protein sequences will be written to this file, in a
#' FASTA format.
#' @param nuc.file If specified, nucleotide version of each protein is written to this file (default
#' \code{NULL}).
#' @param closed.ends Logical, if \code{TRUE} genes are not allowed to run off edges (default \code{TRUE}).
#' @param motif.scan Logical, if \code{TRUE} forces motif scan instead of Shine-Dalgarno trainer (default 
#' \code{FALSE}).
#' 
#' @details This function sets up a call to the software Prodigal (Hyatt et al, 2009). This software is
#' designed to find coding genes in prokaryote genomes. It runs fast and has obtained very good results in
#' tests among the automated gene finders. The options used as default here are believed to be the best for
#' pan-genomic analyses.
#' 
#' @return The call to Prodigal produces a FASTA formatted file with predicted protein sequences, and
#' if \samp{nuc.file} is specified, a similar file with nucleotide sequences. See \code{\link[microseq]{readFasta}}
#' for how to read such files into R.
#' 
#' @references Hyatt, D., Chen, G., LoCascio, P.F., Land, M.L., Larimer, F.W., Hauser, L.J. (2009).
#' Prodigal: prokaryotic gene recognition and translation initiation site identification, BMC Bioinformatics,
#' 11:119.
#' 
#' @note The Prodigal software must be installed on the system for this function to work, i.e. the command 
#' \samp{system("prodigal")} (no version numbers!) must be recognized as a valid command if you run it 
#' in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{entrezDownload}}.
#' 
#' @examples 
#' \dontrun{
#' # Using a small genome file in this package
#' # We need to uncompress it first...
#' extdata.path <- file.path(path.package("micropan"),"extdata")
#' filenames <- "Mpneumoniae_309_genome.fsa"
#' pth <- lapply( file.path( extdata.path, paste( filenames, ".xz", sep="" ) ), xzuncompress )
#' 
#' # Calling Prodigal, and using a similar name (_genome replaced by _protein) in output
#' prodigalPredict( file.path(extdata.path,filenames), gsub("_genome","_protein",filenames) )
#' 
#' # ...and compressing the genome-file again...
#' pth <- lapply( file.path( extdata.path, filenames ), xzcompress )
#' }
#' 
#' @export
prodigalPredict <- function( genome.file, prot.file, nuc.file=NULL, closed.ends=TRUE, motif.scan=FALSE ){
  command <- paste( "prodigal -i ", genome.file, " -a ", prot.file, " -o prodigal.txt -q", sep="" )
  if( !is.null( nuc.file ) ){
    command <- paste( command, " -d ", nuc.file, sep="" )
  }
  if( closed.ends ){
    command <- paste( command, " -c ", sep="" )
  }
  if( motif.scan ){
    command <- paste( command, " -n ", sep="" )
  }
  system( command )
  file.remove( "prodigal.txt" )
  return( paste( "Prodigal predictions in", prot.file ) )
}
