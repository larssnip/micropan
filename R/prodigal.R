#' @name prodigalPredict
#' @title Gene predictions using Prodigal
#' 
#' @description Finds coding genes in a genome using the Prodigal software.
#' 
#' @param genome.file Name of a FASTA formatted file with all the DNA sequences for a genome (chromosomes,
#' plasmids, contigs etc.).
#' @param prot.file If specified, amino acid sequence of each protein is written to this FASTA file.
#' @param nuc.file If specified, nucleotide sequence of each protein is written to this FASTA file.
#' @param closed.ends Logical, if \code{TRUE} genes are not allowed to run off edges (default \code{TRUE}).
#' @param motif.scan Logical, if \code{TRUE} forces motif scan instead of Shine-Dalgarno trainer (default 
#' \code{FALSE}).
#' 
#' @details This function sets up a call to the software Prodigal (Hyatt et al, 2009). This software is
#' designed to find coding genes in prokaryote genomes. It runs fast and has obtained very good results in
#' tests among the automated gene finders. The options used as default here are believed to be the best for
#' pan-genomic analyses.
#' 
#' @return A \code{gff.table} with the metadata for all predicted genes (see \code{\link{readGFF}}). If
#' \code{prot.file} is specified, a FASTA formatted file with predicted protein sequences are also produced. If
#' \code{nuc.file} is specified, a similar file with nucleotide sequences is also produced.
#' 
#' @references Hyatt, D., Chen, G., LoCascio, P.F., Land, M.L., Larimer, F.W., Hauser, L.J. (2010).
#' Prodigal: prokaryotic gene recognition and translation initiation site identification, BMC Bioinformatics,
#' 11:119.
#' 
#' @note The Prodigal software must be installed on the system for this function to work, i.e. the command 
#' \samp{system("prodigal -h")} (no version numbers!) must be recognized as a valid command if you run it 
#' in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{readGFF}}.
#' 
#' @examples 
#' \dontrun{
#' # Using a genome file in this package
#' extdata <- file.path(path.package("micropan"),"extdata")
#' genome.file <- "Mpneumoniae_309_genome.fsa"
#' 
#' # We need to uncompress it first...
#' xzuncompress(file.path(extdata,paste(genome.file,".xz",sep="")))
#'
#' # Calling Prodigal, and writing all predicted proteins to a file as well
#' gff.table <- prodigal(file.path(extdata,genome.file),
#'                       prot.file=gsub("_genome","_protein",genome.file))
#' 
#' # ...and compressing the genome file again...
#' xzcompress(file.path(extdata,genome.file))
#' }
#' 
#' @export prodigal
#' 
prodigal <- function( genome.file, prot.file=NULL, nuc.file=NULL, closed.ends=TRUE, motif.scan=FALSE ){
  cmd <- paste( "prodigal -i ", genome.file, " -f gff -o prodigal.gff -q", sep="" )
  if( !is.null( prot.file ) ){
    cmd <- paste( cmd, " -a ", prot.file, sep="" )
  }
  if( !is.null( nuc.file ) ){
    cmd <- paste( cmd, " -d ", nuc.file, sep="" )
  }
  if( closed.ends ){
    cmd <- paste( cmd, " -c ", sep="" )
  }
  if( motif.scan ){
    cmd <- paste( cmd, " -n ", sep="" )
  }
  system( cmd )
  gff.table <- readGFF( "prodigal.gff" )
  file.remove( "prodigal.gff" )
  return(gff.table)
}
