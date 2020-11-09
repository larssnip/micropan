#' @name extractPanGenes
#' @title Extracting genes of same prevalence
#' 
#' @description Based on a clustering of genes, this function extracts the genes
#' occurring in the same number of genomes.
#' 
#' @param clustering Named vector of clustering 
#' @param N.genomes Vector specifying the number of genomes the genes should be in
#' 
#' @details Pan-genome studies focus on the gene families obtained by some clustering,
#' see \code{\link{bClust}} or \code{\link{dClust}}. This function will extract the individual genes from
#' each genome belonging to gene families found in \code{N.genomes} genomes specified by the user.
#' Only the sequence tag for each gene is extracted, but the sequences can be added easily, see examples
#' below.
#' 
#' @return A table with columns
#' \itemize{
#'   \item cluster. The gene family (integer)
#'   \item seq_tag. The sequence tag identifying each sequence (text)
#'   \item N_genomes. The number of genomes in which it is found (integer)
#' }
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{bClust}}, \code{\link{dClust}}, \code{\link{geneFamilies2fasta}}.
#' 
#' @examples
#' # Loading clustering data in this package
#' data(xmpl.bclst)
#' 
#' # Finding genes in 5 genomes
#' core.tbl <- extractPanGenes(xmpl.bclst, N.genomes = 5)
#' #...or in a single genome
#' orfan.tbl <- extractPanGenes(xmpl.bclst, N.genomes = 1)
#' 
#' \dontrun{
#' # To add the sequences, assume all protein fasta files are in a folder named faa:
#' lapply(list.files("faa", full.names = T), readFasta) %>% 
#'   bind_rows() %>% 
#'   mutate(seq_tag = word(Header, 1, 1)) %>% 
#'   right_join(orfan.tbl, by = "seq_tag") -> orfan.tbl
#'   # The resulting table can be written to fasta file directly using writeFasta()
#'   # See also geneFamilies2fasta()
#' }
#' 
#' @importFrom dplyr distinct group_by summarize filter %>% mutate select arrange
#' @importFrom tibble tibble
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' 
#' @export extractPanGenes
#' 
extractPanGenes <- function(clustering, N.genomes = 1:2){
  tibble(cluster = clustering,
         seq_tag = names(clustering),
         genome_id = str_extract(names(clustering), "GID[0-9]+")) -> tbl
  names(tbl$cluster) <- NULL
  tbl %>% 
    distinct(cluster, genome_id) %>% 
    group_by(cluster) %>% 
    summarize(n.genomes = n()) %>% 
    filter(n.genomes %in% N.genomes) -> trg.tbl
  
  out.tbl <- NULL
  for(i in 1:length(N.genomes)){
    idx <- which(trg.tbl$n.genomes == N.genomes[i])
    tbl %>% 
      filter(cluster %in% trg.tbl$cluster[idx]) %>% 
      mutate(N_genomes = N.genomes[i]) %>% 
      select(cluster, seq_tag, N_genomes) %>% 
      arrange(cluster, seq_tag) %>% 
      bind_rows(out.tbl) -> out.tbl
  }
  return(out.tbl)
}








#' @name geneFamilies2fasta
#' @title Write gene families to files
#' 
#' @description Writes specified gene families to separate fasta files.
#' 
#' @param pangene.tbl A table listing gene families (clusters). 
#' @param fasta.folder The folder containing the fasta files with all sequences.
#' @param out.folder The folder to write to.
#' @param file.ext The file extension to recognize the fasta files in \code{fasta.folder}.
#' @param verbose Logical to allow text ouput during processing
#' 
#' @details The argument \code{pangene.tbl} should be produced by \code{\link{extractPanGenes}} in order to
#' contain the columns \code{cluster}, \code{seq_tag} and \code{N_genomes} required by this function. The
#' files in \code{fasta.folder} must have been prepared by \code{\link{panPrep}} in order to have the proper
#' sequence tag information. They may contain protein sequences or DNA sequences.
#' 
#' If you already added the \code{Header} and \code{Sequence} information to \code{pangene.tbl} these will be
#' used instead of reading the files in \code{fasta.folder}, but a warning is issued.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{extractPanGenes}}, \code{\link{writeFasta}}.
#' 
#' @examples
#' # Loading clustering data in this package
#' data(xmpl.bclst)
#' 
#' # Finding genes in 1,..,5 genomes (all genes)
#' all.tbl <- extractPanGenes(xmpl.bclst, N.genomes = 1:5)
#' 
#' \dontrun{
#' # All protein fasta files are in a folder named faa, and we write to the current folder:
#' clusters2fasta(all.tbl, fasta.folder = "faa", out.folder = ".")
#' 
#' # use pipe, write to folder "orfans"
#' extractPanGenes(xmpl.bclst, N.genomes = 1)) %>% 
#'   geneFamilies2fasta(fasta.folder = "faa", out.folder = "orfans")
#' }
#' 
#' @importFrom dplyr bind_rows distinct filter %>% mutate select right_join
#' @importFrom microseq writeFasta readFasta
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' 
#' @export geneFamilies2fasta
#' 
geneFamilies2fasta <- function(pangene.tbl, fasta.folder, out.folder, file.ext = "fasta$|faa$|fna$|fa$",
                           verbose = TRUE){
  has.Header <- exists("Header", pangene.tbl)
  has.Sequence <- exists("Sequence", pangene.tbl)
  if(has.Header & has.Sequence){
    warning("pangene.tbl already has sequences, ignoring fasta.folder")
  } else {
    if(has.Header) pangene.tbl %>% select(-Header) -> pangene.tbl
    if(has.Sequence) pangene.tbl %>% select(-Sequence) -> pangene.tbl
    fasta.files <- list.files(fasta.folder, pattern = file.ext, full.names = T)
    if(length(fasta.files) == 0) stop("Found no fasta files in fasta.folder")
    if(verbose) cat("clusters2fasta:\n   found", length(fasta.files), "fasta files...\n")
    lapply(fasta.files, readFasta) %>% 
      bind_rows() %>% 
      mutate(seq_tag = word(Header, 1, 1)) %>% 
      right_join(pangene.tbl, by = "seq_tag") -> pangene.tbl
    if(sum(is.na(pangene.tbl$Sequence))) stop("Cannot bind these sequences to the supplied genes, mismatching sequence tags")
  }
  pangene.tbl %>% 
    distinct(cluster, N_genomes) %>% 
    by(1:nrow(.), function(x){str_c("Genome=", x[2], "_Cluster=", x[1], ".fasta")}) %>% 
    as.character() -> filenames
  ucls <- unique(pangene.tbl$cluster)
  if(verbose) cat("   writing", length(ucls), "clusters to files...\n")
  for(i in 1:length(ucls)){
    pangene.tbl %>% 
      filter(cluster == ucls[i]) %>% 
      writeFasta(out.file = file.path(out.folder, filenames[i]))
  }
}