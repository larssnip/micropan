#' @name panPrep
#' @title Preparing FASTA files for pan-genomics
#' 
#' @description Preparing a FASTA file before starting comparisons of sequences.
#' 
#' @usage panPrep(in.file, genome_id, out.file, protein = TRUE, min.length = 10, discard = "")
#' 
#' @param in.file The name of a FASTA formatted file with protein or nucleotide sequences for coding
#' genes in a genome.
#' @param genome_id The Genome Identifier, see below.
#' @param out.file Name of file where the prepared sequences will be written.
#' @param protein Logical, indicating if the \samp{in.file} contains protein (\code{TRUE}) or
#' nucleotide (\code{FALSE}) sequences.
#' @param min.length Minimum sequence length
#' @param discard A text, a regular expression, and sequences having a match against this in their
#' headerline will be discarded.
#' 
#' @details This function will read the \code{in.file} and produce another, slightly modified, FASTA file
#' which is prepared for the comparisons using \code{\link{blastpAllAll}}, \code{\link{hmmerScan}}
#' or any other method.
#' 
#' The main purpose of \code{\link{panPrep}} is to make certain every sequence is labeled with a tag
#' called a \samp{genome_id} identifying the genome from which it comes. This text contains the text
#' \dQuote{GID} followed by an integer. This integer can be any integer as long as it is unique to every
#' genome in the study. If a genome has the text \dQuote{GID12345} as identifier, then the
#' sequences in the file produced by \code{\link{panPrep}} will have headerlines starting with
#' \dQuote{GID12345_seq1}, \dQuote{GID12345_seq2}, \dQuote{GID12345_seq3}...etc. This makes it possible
#' to quickly identify which genome every sequence belongs to.
#' 
#' The \samp{genome_id} is also added to the file name specified in \samp{out.file}. For this reason the
#' \samp{out.file} must have a file extension containing letters only. By convention, we expect FASTA
#' files to have one of the extensions \samp{.fsa}, \samp{.faa}, \samp{.fa} or \samp{.fasta}.
#' 
#' \code{\link{panPrep}} will also remove sequences shorter than \code{min.length}, removing stop codon
#' symbols (\samp{*}), replacing alien characters with \samp{X} and converting all sequences to upper-case.
#' If the input \samp{discard} contains a regular expression, any sequences having a match to this in their
#' headerline are also removed. Example: If we use the \code{prodigal} software (see \code{\link[microseq]{findGenes}})
#' to find proteins in a genome, partially predicted genes will have the text \samp{partial=10} or
#' \samp{partial=01} in their headerline. Using \samp{discard= "partial=01|partial=10"} will remove
#' these from the data set.
#' 
#' @return This function produces a FASTA formatted sequence file, and returns the name of this file.
#' 
#' @author Lars Snipen and Kristian Liland.
#' 
#' @seealso \code{\link{hmmerScan}}, \code{\link{blastpAllAll}}.
#' 
#' @examples
#' # Using a protein file in this package
#' # We need to uncompress it first...
#' pf <- file.path(path.package("micropan"),"extdata","xmpl.faa.xz")
#' prot.file <- tempfile(fileext = ".xz")
#' ok <- file.copy(from = pf, to = prot.file)
#' prot.file <- xzuncompress(prot.file)
#' 
#' # Prepping it, using the genome_id "GID123"
#' prepped.file <- panPrep(prot.file, genome_id = "GID123", out.file = tempfile(fileext = ".faa"))
#' 
#' # Reading the prepped file
#' prepped <- readFasta(prepped.file)
#' head(prepped)
#' 
#' # ...and cleaning...
#' ok <- file.remove(prot.file, prepped.file)
#' 
#' @importFrom microseq readFasta writeFasta
#' @importFrom dplyr mutate filter %>% n
#' @importFrom stringr str_remove_all str_length str_c str_extract str_replace
#' @importFrom rlang .data
#' 
#' @export panPrep
#' 
panPrep <- function(in.file, genome_id, out.file, protein = TRUE, min.length = 10, discard = ""){
  if(protein){
    alien <- "[^ARNDCQEGHILKMFPSTWYV]"
  } else {
    alien <- "[^ACGT]"
  }
  readFasta(normalizePath(in.file)) %>% 
    mutate(Sequence = toupper(.data$Sequence)) %>% 
    mutate(Sequence = str_remove_all(.data$Sequence, "\\*")) %>% 
    mutate(Length = str_length(.data$Sequence)) %>% 
    filter(.data$Length >= min.length) %>% 
    mutate(Sequence = str_replace_all(.data$Sequence, pattern = alien, "X")) %>% 
    mutate(Header = str_c(genome_id, "_seq", 1:n(), " ", .data$Header)) -> fdta
  if(str_length(discard) > 0){
    fdta %>%
      filter(!str_detect(.data$Header, pattern = discard)) -> fdta
  }
  out.file <- normalizePath(out.file)
  fext <- str_extract(out.file, "\\.[a-zA-Z]+$")
  out.file <- str_replace(out.file, fext, str_c("_", genome_id, fext))
  writeFasta(fdta, out.file = out.file)
  return(out.file)
}