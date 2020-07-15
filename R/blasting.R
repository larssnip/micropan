#' @name blastpAllAll
#' @title Making BLAST search all against all genomes
#' 
#' @description Runs a reciprocal all-against-all BLAST search to look for similarity of proteins
#' within and across genomes.
#' 
#' @param prot.files A vector with FASTA filenames.
#' @param out.folder The folder where the result files should end up.
#' @param e.value The chosen E-value threshold in BLAST.
#' @param job An integer to separate multiple jobs.
#' @param start.at An integer to specify where in the file-list to start BLASTing.
#' @param threads The number of CPU's to use.
#' @param verbose Logical, if \code{TRUE} some text output is produced to monitor the progress.
#' 
#' @details A basic step in pangenomics and many other comparative studies is to cluster proteins into
#' groups or families. One commonly used approach is based on BLASTing. This function uses the
#' \samp{blast+} software available for free from NCBI (Camacho et al, 2009). More precisely, the blastp
#' algorithm with the BLOSUM45 scoring matrix and all composition based statistics turned off.
#' 
#' A vector listing FASTA files of protein sequences is given as input in \samp{prot.files}. These files
#' must have the genome_id in the first token of every header, and in their filenames as well, i.e. all input
#' files should first be prepared by \code{\link{panPrep}} to ensure this. Note that only protein sequences
#' are considered here. If your coding genes are stored as DNA, please translate them to protein prior to
#' using this function, see \code{\link[microseq]{translate}}.
#' 
#' In the first version of this package we used reciprocal BLASTing, i.e. we computed both genome A against
#' B and B against A. This may sometimes produce slightly different results, but in reality this is too
#' costly compared to its gain, and we now only make one of the above searches. This basically halves the
#' number of searches. This step is still very time consuming for larger number of genomes. Note that the 
#' protein files are sorted by the genome_id (part of filename) inside this function. This is to ensure a 
#' consistent ordering irrespective of how they are enterred.
#' 
#' For every pair of genomes a result file is produced. If two genomes have genome_id's \samp{GID111},
#' and \samp{GID222} then the result file \samp{GID222_vs_GID111.txt} will
#' be found in \samp{out.folder} after the completion of this search. The last of the two genome_id is always
#' the first in alphabetical order of the two.
#' 
#' The \samp{out.folder} is scanned for already existing result files, and \code{\link{blastpAllAll}} never
#' overwrites an existing result file. If a file with the name \samp{GID111_vs_GID222.txt} already exists in
#' the \samp{out.folder}, this particular search is skipped. This makes it possible to run multiple jobs in
#' parallell, writing to the same \samp{out.folder}. It also makes it possible to add new genomes, and only
#' BLAST the new combinations without repeating previous comparisons. 
#' 
#' This search can be slow if the genomes contain many proteins and it scales quadratically in the number of
#' input files. It is best suited for the study of a smaller number of genomes. By
#' starting multiple R sessions, you can speed up the search by running \code{\link{blastpAllAll}} from each R
#' session, using the same \samp{out.folder} but different integers for the \code{job} option. At the same
#' time you may also want to start the BLASTing at different places in the file-list, by giving larger values
#' to the argument \code{start.at}. This is 1 by default, i.e. the BLASTing starts at the first protein file.
#' If you are using a multicore computer you can also increase the number of CPUs by increasing \code{threads}.
#' 
#' The result files are tab-separated text files, and can be read into R, but more
#' commonly they are used as input to \code{\link{bDist}} to compute distances between sequences for subsequent
#' clustering.
#' 
#' @return The function produces a result file for each pair of files listed in \samp{prot.files}.
#' These result files are located in \code{out.folder}. Existing files are never overwritten by
#' \code{\link{blastpAllAll}}, if you want to re-compute something, delete the corresponding result files first.
#' 
#' @references Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L.
#' (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10:421.
#' 
#' @note The \samp{blast+} software must be installed on the system for this function to work, i.e. the command
#' \samp{system("makeblastdb -help")} must be recognized as valid commands if you
#' run them in the Console window.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panPrep}}, \code{\link{bDist}}.
#' 
#' @examples 
#' \dontrun{
#' # This example requires the external BLAST+ software
#' # Using protein files in this package
#' pf <- file.path(path.package("micropan"), "extdata",
#'                 str_c("xmpl_GID", 1:3, ".faa.xz"))
#' 
#' # We need to uncompress them first...
#' prot.files <- tempfile(fileext = c("_GID1.faa.xz","_GID2.faa.xz","_GID3.faa.xz"))
#' ok <- file.copy(from = pf, to = prot.files)
#' prot.files <- unlist(lapply(prot.files, xzuncompress))
#' 
#' # Blasting all versus all
#' out.dir <- "."
#' blastpAllAll(prot.files, out.folder = out.dir)
#' 
#' # Reading results, and computing blast.distances
#' blast.files <- list.files(out.dir, pattern = "GID[0-9]+_vs_GID[0-9]+.txt")
#' blast.distances <- bDist(file.path(out.dir, blast.files))
#' 
#' # ...and cleaning...
#' ok <- file.remove(prot.files)
#' ok <- file.remove(file.path(out.dir, blast.files))
#' }
#' 
#' @importFrom stringr str_extract str_c
#' 
#' @export blastpAllAll
blastpAllAll <- function(prot.files, out.folder, e.value = 1, job = 1, threads = 1, start.at = 1, verbose = TRUE){
  if(available.external("blast+")){
    N <- length(prot.files)
    genome_id <- str_extract(prot.files, "GID[0-9]+")
    prot.files <- prot.files[order(genome_id)]
    genome_id <- genome_id[order(genome_id)]
    prot.files <- prot.files[start.at:N]
    genome_id <- genome_id[start.at:N]
    N <- length(prot.files)
    if(.Platform$OS.type == "windows"){
      outfmt <- "\"6 qseqid sseqid evalue bitscore\""
    } else {
      outfmt <- "'6 qseqid sseqid evalue bitscore'"
    }
    out.folder <- normalizePath(out.folder)
    file.tbl <- data.frame(Dbase = rep("", (N^2 + N)/2),
                           Query = rep("", (N^2 + N)/2),
                           Res.file = rep("", (N^2 + N)/2),
                           stringsAsFactors = FALSE)
    cc <- 1
    for(i in 1:N){
      for(j in i:N){
        file.tbl$Dbase[cc] <- prot.files[i]
        file.tbl$Query[cc] <- prot.files[j]
        file.tbl$Res.file[cc] <- str_c(genome_id[j], "_vs_", genome_id[i], ".txt")
        cc <- cc + 1
      }
    }
    existing.files <- list.files(out.folder, pattern = "txt$")
    file.tbl %>% 
      filter(!(.data$Res.file %in% existing.files)) -> file.tbl
    if(nrow(file.tbl) > 0){
      dbases <- unique(file.tbl$Dbase)
      for(i in 1:length(dbases)){
        log.fil <- file.path(out.folder, str_c("log", job, ".txt"))
        db.fil <- file.path(out.folder, str_c("blastDB", job))
        if(verbose) cat("blastpAllAll: Making BLAST database of", dbases[i], "\n")
        system(paste("makeblastdb -logfile",log.fil, "-dbtype prot -out", db.fil, "-in", dbases[i]))
        file.tbl %>% 
          filter(.data$Dbase == dbases[i]) -> tbl
        for(j in 1:nrow(tbl)){
          out.file <- file.path(out.folder, tbl$Res.file[j])
          if(!file.exists(out.file)){
            if(verbose) cat("   ", tbl$Res.file[j], "\n")
            cmd <- paste("blastp",
                         "-matrix BLOSUM45",
                         "-evalue", e.value,
                         "-num_threads", threads,
                         "-comp_based_stats", "F",
                         "-num_alignments", 1000,
                         "-outfmt", outfmt,
                         "-query", tbl$Query[j],
                         "-db", db.fil,
                         "-out", out.file)
            system(cmd)
          }
        }
      }
      ok <- file.remove(list.files(out.folder, pattern = str_c("blastDB", job), full.names = T))
      ok <- file.remove(log.fil)
      if(file.exists(str_c(log.fil, ".perf"))) ok <- file.remove(str_c(log.fil, ".perf"))
    }
    invisible(TRUE)
  }
}
