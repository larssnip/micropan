#' @name panPrep
#' @title Preparing FASTA files for pan-genomics
#' 
#' @description Preparing a FASTA file before starting comparisons of sequences in a pan-genome study.
#' 
#' @param in.file The name of a FASTA formatted file with protein or nucleotide sequences for coding
#' genes in a genome.
#' @param GID.tag The Genome IDentifier tag, see below.
#' @param out.file Name of file where the prepared sequences will be written.
#' @param protein Logical, indicating if the \samp{in.file} contains protein (\code{TRUE}) or
#' nucleotide (\code{FALSE}) sequences.
#' @param discard A text, a regular expression, and sequences having a match against this in their
#' header text will be discarded.
#' 
#' @details This function will read a FASTA file and produce another, slightly modified, FASTA file
#' which is prepared for genome-wise comparisons using \code{\link{blastAllAll}}, \code{\link{hmmerScan}}
#' or any other method. 
#' 
#' The main purpose of \code{\link{panPrep}} is to make certain every sequence is labeled with a tag
#' called a \samp{GID.tag} (Genome IDentifier tag) identifying the genome. This text contains the text
#' \dQuote{GID} followed by an integer. This integer can be any integer as long as it is unique to every
#' genome in the study. It can typically be the BioProject number or any other integer that is uniquely
#' related to a specific genome. If a genome has the text \dQuote{GID12345} as identifier, then the
#' sequences in the file produced by \code{\link{panPrep}} will have headerlines starting with
#' \dQuote{GID12345_seq1}, \dQuote{GID12345_seq2}, \dQuote{GID12345_seq3}...etc. This makes it possible
#' to quickly identify which genome every sequence belongs to.
#' 
#' The \samp{GID.tag} is also added to the file name specified in \samp{out.file}. For this reason the
#' \samp{out.file} must have a file extension containing letters only. By convention, we expect FASTA
#' files to have one of the extensions \samp{.fsa}, \samp{.faa}, \samp{.fa} or \samp{.fasta}.
#' 
#' \code{\link{panPrep}} will also remove very short sequences (< 10 amino acids), removing stop codon
#' symbols (\samp{*}), replacing alien characters with \samp{X} and converting all sequences to upper-case.
#' If the input \samp{discard} contains a regular expression, any sequences having a match to this in their
#' headerline are also removed. Example: If we use \code{\link{prodigalPredict}} to find proteins in a
#' genome, partially predicted genes will have the text \samp{partial=10} or \samp{partial=01} in their
#' headerline. Using \samp{discard="partial=01|partial=10"} will remove these from the data set.
#' 
#' @return This function produces a FASTA formatted sequence file.
#' 
#' @author Lars Snipen and Kristian Liland.
#' 
#' @seealso \code{\link{hmmerScan}}, \code{\link{blastAllAll}}.
#' 
#' @examples
#' \dontrun{
#' # Using a protein file in the micropan package
#' extdata <- file.path(path.package("micropan"),"extdata")
#' prot.file <- "Mpneumoniae_309_protein.fsa"
#' 
#' # We need to uncompress it first...
#' xzuncompress(file.path(extdata,paste(prot.file,".xz",sep="")))
#' 
#' # Prepping it, using the GID.tag "GID123"
#' panPrep(file.path(extdata,prot.file),GID.tag="GID123","Mpneumoniae_309.fsa") 
#' # ...should produce a FASTA file named Mpneumoniae_309_GID123.fsa
#' 
#' # ...and compress the input file again...
#' xzcompress(file.path(extdata,prot.file))
#' }
#' 
#' @importFrom microseq readFasta writeFasta
#' 
#' @export panPrep
#' 
panPrep <- function( in.file, GID.tag, out.file, protein=TRUE, discard=NA ){
  if( protein ){
    alien <- "[^ARNDCQEGHILKMFPSTWYV]"
    minl <- 10
  } else {
    alien <- "[^ACGT]"
    minl <- 30
  }
  fdta <- readFasta( in.file )
  if( !is.na( discard ) ) fdta <- fdta[grep( discard, fdta$Header, invert=T ),]
  fdta$Sequence <- toupper( gsub( "\\*", "", fdta$Sequence ) )
  fdta <- fdta[which( nchar( fdta$Sequence ) >= minl ),]
  fdta$Sequence <- gsub( alien, "X", fdta$Sequence )
  nseq <- dim( fdta )[1]
  tok1 <- paste( GID.tag, paste( "seq", 1:nseq, sep="" ), sep="_" )
  fdta$Header <- paste( tok1, fdta$Header )
  fext <- unlist( microseq::gregexpr( "\\.[a-zA-Z]+$", out.file, extract=T ) )
  fname <- paste( gsub( "\\.[a-zA-Z]+$", "", out.file ), "_", GID.tag, fext, sep="" )
  writeFasta( fdta, out.file=fname )
}