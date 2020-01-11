#' @name bDist
#' @title Computes distances between sequences based on BLAST results
#' 
#' @description Reads a complete set of result files from a BLAST search and
#' computes distance between all sequences based on the BLAST bit-score.
#' 
#' @param blast.files A text vector of filenames.
#' @param e.value A threshold E-value to immediately discard (very) poor BLAST alignments.
#' @param verbose Logical, indicating if textual output should be given to monitor the progress.
#' 
#' @details Each input file must be a BLAST result file produced by \code{\link{blastpAllAll}}.
#' 
#' Setting a small \samp{e.value} threshold can speed up the computation and reduce the size of the
#' output, but you may loose some alignments that could produce smallish distances for short sequences.
#' 
#' The distance computed is based on alignment bitscores. Assume the alignment of query A against hit B
#' has a bitscore of S(A,B). Then, we also typically find an alignment of query B against hit A with
#' bitscore S(B,A). The maximum of these two are used as the bitscore for this pair, and in virtually
#' all cases S(A,B)=S(B,A) anyway. The distance is D(A,B)=1-2*S(A,B)/(S(A,A)+S(B,B)) where S(A,A) (or S(B,B))
#' is the bitscore of aligning A (or B) against itself. A distance of
#' 0.0 means A and B are identical. The maximum possible distance is 1.0, meaning there is no BLAST hit
#' found either way.
#'
#' This distance should not be interpreted as lack of identity! A distance of 0.0 means 100\% identity,
#' but a distance of 0.25 does \emph{not} mean 75\% identity. It has some resemblance to an evolutinary
#' (raw) distance, but since it is based on protein alignments, the type of mutations plays a significant
#' role, not only the number of mutations.
#' 
#' @return The function returns a \samp{tibble} with columns \samp{Query}, \samp{Hit}, \samp{Bitscore}
#' and \samp{Distance}. Each row corresponds to a pair of sequences having at least one BLAST hit between
#' them. All pairs \emph{not} listed in the output have distance 1.0 between them.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{blastpAllAll}}, \code{\link{bClust}}, \code{\link{isOrtholog}}.
#' 
#' @examples
#' # Using BLAST result files in this package...
#' prefix <- c("GID1_vs_GID1.txt",
#'             "GID1_vs_GID2.txt",
#'             "GID1_vs_GID3.txt",
#'             "GID2_vs_GID2.txt",
#'             "GID2_vs_GID3.txt",
#'             "GID3_vs_GID3.txt")
#' bf <- file.path(path.package("micropan"), "extdata", str_c(prefix, ".xz"))
#' 
#' # We need to uncompress them first...
#' blast.files <- tempfile(pattern = prefix, fileext = ".xz")
#' ok <- file.copy(from = bf, to = blast.files)
#' blast.files <- unlist(lapply(blast.files, xzuncompress))
#' 
#' # Computing pairwise distances
#' blast.dist <- bDist(blast.files)
#' 
#' # ...and cleaning...
#' ok <- file.remove(blast.files)
#' 
#' # See also example for blastpAllAll
#' 
#' @importFrom tibble tibble
#' @importFrom stringr str_extract
#' @importFrom readr read_delim
#' 
#' @export
bDist <- function(blast.files, e.value = 1, verbose = TRUE){
  if(verbose) cat("bDist:\n")
  gids <- str_extract_all(blast.files, "GID[0-9]+", simplify = T)
  self.idx <- which(gids[,1] == gids[,2])
  if(verbose) cat("...reading", length(self.idx), "self alignments...\n")
  
  slf.tbl <- NULL
  max.tbl <- NULL
  for(i in 1:length(self.idx)){
    suppressMessages(read_delim(blast.files[self.idx[i]], delim = "\t", trim_ws = T,
                               col_names = c("Query", "Hit", "Evalue", "Bitscore"))) %>% 
      filter(Evalue <= e.value) %>% 
      arrange(desc(Bitscore)) %>% 
      mutate(Pair = sortPaste(Query, Hit)) %>% 
      distinct(Pair, .keep_all = T) %>% 
      select(-Evalue, -Pair) -> tbl
    tbl %>%
      filter(Query == Hit) %>% 
      bind_rows(max.tbl) -> max.tbl
    idx.q <- match(tbl$Query, max.tbl$Query)
    idx.h <- match(tbl$Hit,   max.tbl$Query)
    tbl %>%
      mutate(Distance = 1 - (2 * Bitscore) / (max.tbl$Bitscore[idx.q] + max.tbl$Bitscore[idx.h])) %>% 
      bind_rows(slf.tbl) -> slf.tbl
  }
  if(verbose) cat("...found BLAST results for", nrow(slf.tbl), "unique sequences...\n")

  if(verbose) cat("...reading remaining alignments...\n")
  crss.tbl <- NULL
  for(i in 1:length(blast.files)){
    if(!(i %in% self.idx)){
      suppressMessages(read_delim(blast.files[i], delim = "\t", trim_ws = T,
                                  col_names = c("Query", "Hit", "Evalue", "Bitscore"))) %>% 
        filter(Evalue <= e.value) %>% 
        arrange(desc(Bitscore)) %>% 
        mutate(Pair = sortPaste(Query, Hit)) %>% 
        distinct(Pair, .keep_all = T) %>% 
        select(-Evalue, -Pair) -> tbl
      idx.q <- match(tbl$Query, max.tbl$Query)
      idx.h <- match(tbl$Hit,   max.tbl$Query)
      if(sum(is.na(idx.q)) > 0) stop("Self-alignment lacking for Query", which(is.na(idx.q)), "in blast.file", blast.files[i], "\n")
      if(sum(is.na(idx.h)) > 0) stop("Self-alignment lacking for Hit", which(is.na(idx.h)), "in blast.file", blast.files[i], "\n")
      tbl %>% 
        mutate(Distance = 1 - (2 * Bitscore)/(max.tbl$Bitscore[idx.q] + max.tbl$Bitscore[idx.h])) %>% 
        bind_rows(crss.tbl) -> crss.tbl
    }
  }
  bind_rows(slf.tbl, crss.tbl) %>% 
    arrange(Query, Hit) -> dist.tbl
  return(dist.tbl)
}


# Local function
sortPaste <- function(q, h){
  M <- matrix(c(q, h), ncol = 2, byrow = F)
  pp <- apply(M, 1, function(x){paste(sort(x), collapse = ":")})
  return(pp)
}
