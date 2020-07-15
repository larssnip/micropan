#' @name bDist
#' @title Computes distances between sequences
#' 
#' @description Computes distance between all sequences based on the BLAST bit-scores.
#' 
#' @param blast.files A text vector of BLAST result filenames.
#' @param blast.tbl A table with BLAST results.
#' @param e.value A threshold E-value to immediately discard (very) poor BLAST alignments.
#' @param verbose Logical, indicating if textual output should be given to monitor the progress.
#' 
#' @details The essential input is either a vector of BLAST result filenames (\code{blast.files}) or a
#' table of the BLAST results (\code{blast.tbl}). It is no point in providing both, then \code{blast.tbl} is ignored.
#' 
#' For normal sized data sets (e.g. less than 100 genomes), you would provide the BLAST filenames as the argument
#' \code{blast.files} to this function.
#' Then results are read, and distances are computed. Only if you have huge data sets, you may find it more efficient to 
#' read the files using \code{\link{readBlastSelf}} and \code{\link{readBlastPair}} separately, and then provide as the
#' argument \code{blast.tbl]} the table you get from binding these results. In all cases, the BLAST result files must
#' have been produced by \code{\link{blastpAllAll}}.
#' 
#' Setting a small \samp{e.value} threshold can speed up the computation and reduce the size of the
#' output, but you may loose some alignments that could produce smallish distances for short sequences.
#' 
#' The distance computed is based on alignment bitscores. Assume the alignment of query A against hit B
#' has a bitscore of S(A,B). The distance is D(A,B)=1-2*S(A,B)/(S(A,A)+S(B,B)) where S(A,A) and S(B,B) are
#' the self-alignment bitscores, i.e. the scores of aligning against itself. A distance of
#' 0.0 means A and B are identical. The maximum possible distance is 1.0, meaning there is no BLAST between A and B.
#'
#' This distance should not be interpreted as lack of identity! A distance of 0.0 means 100\% identity,
#' but a distance of 0.25 does \emph{not} mean 75\% identity. It has some resemblance to an evolutinary
#' (raw) distance, but since it is based on protein alignments, the type of mutations plays a significant
#' role, not only the number of mutations.
#' 
#' @return The function returns a table with columns \samp{Dbase}, \samp{Query}, \samp{Bitscore}
#' and \samp{Distance}. Each row corresponds to a pair of sequences (Dbase and Query sequences) having at least
#' one BLAST hit between
#' them. All pairs \emph{not} listed in the output have distance 1.0 between them.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{blastpAllAll}}, \code{\link{readBlastSelf}}, \code{\link{readBlastPair}},
#' \code{\link{bClust}}, \code{\link{isOrtholog}}.
#' 
#' @examples
#' # Using BLAST result files in this package...
#' prefix <- c("GID1_vs_GID1_",
#'             "GID2_vs_GID1_",
#'             "GID3_vs_GID1_",
#'             "GID2_vs_GID2_",
#'             "GID3_vs_GID2_",
#'             "GID3_vs_GID3_")
#' bf <- file.path(path.package("micropan"), "extdata", str_c(prefix, ".txt.xz"))
#' 
#' # We need to uncompress them first...
#' blast.files <- tempfile(pattern = prefix, fileext = ".txt.xz")
#' ok <- file.copy(from = bf, to = blast.files)
#' blast.files <- unlist(lapply(blast.files, xzuncompress))
#' 
#' # Computing pairwise distances
#' blast.dist <- bDist(blast.files)
#' 
#' # Read files separately, then use bDist
#' self.tbl <- readBlastSelf(blast.files)
#' pair.tbl <- readBlastPair(blast.files)
#' blast.dist <- bDist(blast.tbl = bind_rows(self.tbl, pair.tbl))
#' 
#' # ...and cleaning...
#' ok <- file.remove(blast.files)
#' 
#' # See also example for blastpAl
#' 
#' @importFrom tibble tibble
#' @importFrom stringr str_extract_all
#' @importFrom dplyr %>% rename filter arrange mutate distinct select bind_rows desc
#' @importFrom utils read.table
#' @importFrom rlang .data
#' 
#' @export bDist
#' 
bDist <- function(blast.files = NULL, blast.tbl = NULL, e.value = 1, verbose = TRUE){
  if(!is.null(blast.files)){
    readBlastSelf(blast.files, e.value = e.value, verbose = verbose) %>% 
      filter(.data$Evalue <= e.value) %>% 
      arrange(desc(.data$Bitscore)) %>% 
      mutate(Pair = sortPaste(.data$Dbase, .data$Query)) %>% 
      distinct(.data$Pair, .keep_all = TRUE) %>% 
      select(.data$Dbase, .data$Query, .data$Bitscore) -> self.tbl
    readBlastPair(blast.files, e.value = e.value, verbose = verbose) %>% 
      filter(.data$Evalue <= e.value) %>% 
      select(.data$Dbase, .data$Query, .data$Bitscore) %>% 
      arrange(desc(.data$Bitscore)) %>% 
      distinct(.data$Dbase, .data$Query, .keep_all = TRUE) %>% 
      bind_rows(self.tbl) -> blast.tbl
  } else if(is.null(blast.tbl)){
    stop("Needs either blast.files or blast.tbl as input")
  } else {
    blast.tbl %>% 
      filter(.data$Evalue <= e.value) %>% 
      select(-.data$Evalue) %>% 
      arrange(desc(.data$Bitscore)) %>% 
      mutate(GIDd = str_extract(.data$Dbase, "GID[0-9]+")) %>% 
      mutate(GIDq = str_extract(.data$Query, "GID[0-9]+")) -> blast.tbl
    blast.tbl %>% 
      filter(.data$GIDd == .data$GIDq) %>% 
      mutate(Pair = sortPaste(.data$Dbase, .data$Query)) %>% 
      distinct(.data$Pair, .keep_all = TRUE) %>% 
      select(-.data$Pair) -> self.tbl
    blast.tbl %>% 
      filter(.data$GIDd != .data$GIDq) %>% 
      bind_rows(self.tbl) %>% 
      select(-.data$GIDd, -.data$GIDq) %>% 
      distinct(.data$Dbase, .data$Query, .keep_all = TRUE) -> blast.tbl
  }
  if(verbose) cat("bDist:\n   ...found", nrow(blast.tbl), "alignments...\n")
  blast.tbl %>% 
    filter(.data$Dbase == .data$Query) -> self.tbl
  if(verbose) cat("   ...where", nrow(self.tbl), "are self-alignments...\n")
  idx.d <- match(blast.tbl$Dbase, self.tbl$Dbase)
  idd <- which(is.na(idx.d))
  if(length(idd) > 0) stop("No self-alignment for sequences: ", str_c(unique(blast.tbl$Dbase[idd]), collapse = ","))
  idx.q <- match(blast.tbl$Query, self.tbl$Query)
  idd <- which(is.na(idx.q))
  if(length(idd) > 0) stop("No self-alignment for sequences: ", str_c(unique(blast.tbl$Query[idd]), collapse = ","))
  
  blast.tbl %>%
    mutate(Distance = 1 - (2 * .data$Bitscore) / (self.tbl$Bitscore[idx.d] + self.tbl$Bitscore[idx.q])) %>% 
    arrange(.data$Dbase, .data$Query) -> dist.tbl
  return(dist.tbl)
}


# Local function
sortPaste <- function(q, h){
  M <- matrix(c(q, h), ncol = 2, byrow = F)
  pp <- apply(M, 1, function(x){paste(sort(x), collapse = ":")})
  return(pp)
}



#' @name readBlastSelf
#' @aliases readBlastSelf readBlastPair
#' @title Reads BLAST result files
#' 
#' @description Reads files from a search with blastpAllAll
#' 
#' @param blast.files A text vector of filenames.
#' @param e.value A threshold E-value to immediately discard (very) poor BLAST alignments.
#' @param verbose Logical, indicating if textual output should be given to monitor the progress.
#' 
#' @details The filenames given as input must refer to BLAST result files produced by \code{\link{blastpAllAll}}.
#' 
#' With \code{readBlastSelf} you only read the self-alignment results, i.e. blasting a genome against itself. With
#' \code{readBlastPair} you read all the other files, i.e. different genomes compared. You may use all blast file
#' names as input to both, they will select the proper files based on their names, e.g. GID1_vs_GID1.txt is read
#' by \code{readBlastSelf} while GID2_vs_GID1.txt is read by \code{readBlastPair}.
#' 
#' Setting a small \samp{e.value} threshold will filter the alignment, and may speed up this and later processing,
#' but you may also loose some important alignments for short sequences.
#' 
#' Both these functions are used by \code{\link{bDist}}. The reason we provide them separately is to allow the user
#' to complete this file reading before calling \code{\link{bDist}}. If you have a huge number of files, a
#' skilled user may utilize parallell processing to speed up the reading. For normal size data sets (e.g. less than 100 genomes)
#' you should probably use \code{\link{bDist}} directly.
#' 
#' @return The functions returns a table with columns \samp{Dbase}, \samp{Query}, \samp{Bitscore}
#' and \samp{Distance}. Each row corresponds to a pair of sequences (a Dbase and a Query sequence) having at least
#' one BLAST hit between
#' them. All pairs \emph{not} listed have distance 1.0 between them. You should normally bind the output from 
#' \code{readBlastSelf} to the ouptut from \code{readBlastPair} and use the result as input to \code{\link{bDist}}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{bDist}}, \code{\link{blastpAllAll}}.
#' 
#' @examples
#' # Using BLAST result files in this package...
#' prefix <- c("GID1_vs_GID1_",
#'             "GID2_vs_GID1_",
#'             "GID3_vs_GID1_",
#'             "GID2_vs_GID2_",
#'             "GID3_vs_GID2_",
#'             "GID3_vs_GID3_")
#' bf <- file.path(path.package("micropan"), "extdata", str_c(prefix, ".txt.xz"))
#' 
#' # We need to uncompress them first...
#' blast.files <- tempfile(pattern = prefix, fileext = ".txt.xz")
#' ok <- file.copy(from = bf, to = blast.files)
#' blast.files <- unlist(lapply(blast.files, xzuncompress))
#' 
#' # Reading self-alignment files, then the other files
#' self.tbl <- readBlastSelf(blast.files)
#' pair.tbl <- readBlastPair(blast.files)
#' 
#' # ...and cleaning...
#' ok <- file.remove(blast.files)
#' 
#' # See also examples for bDist
#' 
#' @importFrom stringr str_extract_all
#' @importFrom dplyr %>% rename filter bind_rows
#' @importFrom utils read.table
#' @importFrom rlang .data
#' 
#' @export readBlastSelf readBlastPair
#' 
readBlastSelf <- function(blast.files, e.value = 1, verbose = TRUE){
  blast.files <- normalizePath(blast.files)
  if(verbose) cat("readBlastSelf:\n   ...received", length(blast.files), "blast-files...\n")
  gids <- str_extract_all(blast.files, "GID[0-9]+", simplify = TRUE)
  self.idx <- which(gids[,1] == gids[,2])
  if(verbose) cat("   ...found", length(self.idx), "self-alignment files...\n")
  lapply(blast.files[self.idx], read.table, header = FALSE, sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE) %>% 
    bind_rows() %>% 
    rename(Dbase = .data$V1, Query = .data$V2, Evalue = .data$V3, Bitscore = .data$V4) %>% 
    filter(.data$Evalue <= e.value) -> self.tbl
  if(verbose) cat("   ...returns", nrow(self.tbl), "alignment results\n")
  return(self.tbl)
}
readBlastPair <- function(blast.files, e.value = 1, verbose = TRUE){
  blast.files <- normalizePath(blast.files)
  if(verbose) cat("readBlastPairs:\n   ...received", length(blast.files), "blast-files...\n")
  gids <- str_extract_all(blast.files, "GID[0-9]+", simplify = TRUE)
  pair.idx <- which(gids[,1] != gids[,2])
  if(verbose) cat("   ...found", length(pair.idx), "alignment files who are NOT self-alignments...\n")
  lapply(blast.files[pair.idx], read.table, header = FALSE, sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE) %>% 
    bind_rows() %>% 
    rename(Dbase = .data$V1, Query = .data$V2, Evalue = .data$V3, Bitscore = .data$V4) %>% 
    filter(.data$Evalue <= e.value) -> pair.tbl
  if(verbose) cat("   ...returns", nrow(pair.tbl), "alignment results\n")
  return(pair.tbl)
}

