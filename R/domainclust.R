#' @name dClust
#' @title Clustering sequences based on domain sequence
#' 
#' @description Proteins are clustered by their sequence of protein domains. A domain sequence is the
#' ordered sequence of domains in the protein. All proteins having identical domain sequence are assigned
#' to the same cluster.
#' 
#' @param hmmer.tbl A \code{tibble} of results from a \code{\link{hmmerScan}} against a domain database.
#' 
#' @details A domain sequence is simply the ordered list of domains occurring in a protein. Not all proteins
#' contain known domains, but those who do will have from one to several domains, and these can be ordered
#' forming a sequence. Since domains can be more or less conserved, two proteins can be quite different in
#' their amino acid sequence, and still share the same domains. Describing, and grouping, proteins by their
#' domain sequence was proposed by Snipen & Ussery (2012) as an alternative to clusters based on pairwise
#' alignments, see \code{\link{bClust}}. Domain sequence clusters are less influenced by gene prediction errors.
#' 
#' The input is a \code{tibble} of the type produced by \code{\link{readHmmer}}. Typically, it is the
#' result of scanning proteins (using \code{\link{hmmerScan}}) against Pfam-A or any other HMMER3 database
#' of protein domains. It is highly reccomended that you remove overlapping hits in \samp{hmmer.tbl} before
#' you pass it as input to \code{\link{dClust}}. Use the function \code{\link{hmmerCleanOverlap}} for this.
#' Overlapping hits are in some cases real hits, but often the poorest of them are artifacts.
#' 
#' @return The output is a numeric vector with one element for each unique sequence in the \samp{Query}
#' column of the input \samp{hmmer.tbl}. Sequences with identical number belong to the same cluster. The
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
#' # HMMER3 result files in this package
#' hf <- file.path(path.package("micropan"), "extdata", 
#'                 str_c("GID", 1:3, "_vs_Pfam-A.hmm.txt.xz"))
#' 
#' # We need to uncompress them first...
#' hmm.files <- tempfile(fileext = rep(".xz", length(hf)))
#' ok <- file.copy(from = hf, to = hmm.files)
#' hmm.files <- unlist(lapply(hmm.files, xzuncompress))
#' 
#' # Reading the HMMER3 results, cleaning overlaps...
#' hmmer.tbl <- NULL
#' for(i in 1:3){
#'   readHmmer(hmm.files[i]) %>% 
#'     hmmerCleanOverlap() %>% 
#'     bind_rows(hmmer.tbl) -> hmmer.tbl
#' }
#' 
#' # The clustering
#' clustering.domains <- dClust(hmmer.tbl)
#' 
#' # ...and cleaning...
#' ok <- file.remove(hmm.files)
#' 
#' @export dClust
#' 
dClust <- function(hmmer.tbl){
  hmmer.tbl %>% 
    arrange(Start) %>% 
    group_by(Query) %>% 
    summarize(Dom.seq = str_c(Hit, collapse = ",")) %>% 
    mutate(Cluster = as.integer(factor(Dom.seq, levels = unique(Dom.seq)))) -> tbl

  dsc <- tbl$Cluster
  names(dsc) <- tbl$Query
  attr(dsc, "cluster.info") <- unique(tbl$Dom.seq)
  return(dsc)
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
#' The input is a \code{tibble} of the type produced by \code{\link{readHmmer}}.
#' 
#' @return A \code{tibble} which is a subset of the input, where some rows may have been deleted to
#' avoid overlapping hits.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{hmmerScan}}, \code{\link{readHmmer}}, \code{\link{dClust}}.
#' 
#' @examples # See the example in the Help-file for dClust.
#' 
#' @importFrom dplyr filter select %>% 
#' 
#' @export hmmerCleanOverlap
#' 
hmmerCleanOverlap <- function(hmmer.tbl){
  qt <- table(hmmer.tbl$Query)
  if(max(qt) > 1){
    multi <- names(qt[qt > 1])
    hmmer.tbl$Keep <- TRUE
    for(i in 1:length(multi)){
      idx <- which(hmmer.tbl$Query == multi[i])
      hmmer.tbl$Keep[idx] <- keeper(hmmer.tbl[idx,])
    }
    hmmer.tbl %>% 
      filter(Keep) %>% 
      select(-Keep) -> hmmer.tbl
  }
  return(hmmer.tbl)
} 



# Local functions
keeper<- function(hmmer.tbl){
  hmmer.tbl$Overlaps <- overlapper(hmmer.tbl)
  while((sum(hmmer.tbl$Keep) > 1) & (sum(hmmer.tbl$Overlaps) > 0)){
    idx <- which(hmmer.tbl$Overlaps)
    idd <- which(hmmer.tbl$Evalue[idx] == max(hmmer.tbl$Evalue[idx]))
    hmmer.tbl$Keep[idx[idd[1]]] <- F
    hmmer.tbl$Overlaps <- overlapper(hmmer.tbl)
  }
  return(hmmer.tbl$Keep)
}
overlapper <- function(hmmer.tbl){
  olaps <- rep(FALSE, nrow(hmmer.tbl))
  idx <- which(hmmer.tbl$Keep)
  if(length(idx) > 1){
    ht <- slice(hmmer.tbl, idx)
    for(i in 1:nrow(ht)){
      ovr <- ((ht$Start[i] <= ht$Stop[-i]) & (ht$Start[i] >= ht$Start[-i])) |
        ((ht$Stop[i] <= ht$Stop[-i]) & (ht$Stop[i] >= ht$Start[-i])) |
        ((ht$Start[i] <= ht$Start[-i]) & (ht$Stop[i] >= ht$Stop[-i]))
      olaps[idx[i]] <- ovr && ovr
    }
  }
  return(olaps)
}
