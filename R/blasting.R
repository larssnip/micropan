#' @name blastpAllAll
#' @title Making BLAST search all against all genomes
#' 
#' @description Runs a reciprocal all-against-all BLAST search to look for similarity of proteins
#' within and across genomes.
#' 
#' @param prot.files A vector with the names of the FASTA files where the protein sequences of
#' each genome is found.
#' @param out.folder The name of the folder where the result files should end up.
#' @param e.value The chosen E-value threshold in BLAST. Default is \samp{e.value = 1}, a smaller value
#' will speed up the search at the cost of less sensitivity.
#' @param job An integer to separate multiple jobs. You may want to run several jobs in parallell,
#' and each job should have different number here to avoid confusion on databases.
#' @param threads The number of CPU's to use.
#' @param verbose Logical, if \code{TRUE} some text output is produced to monitor the progress.
#' 
#' @details A basic step in pangenomics and many other comparative studies is to cluster proteins into
#' groups or families. One commonly used approach is based on BLASTing. This function uses the
#' \samp{blast+} software available for free from NCBI (Camacho et al, 2009). 
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
#' number of searches. This step is still very time consuming for larger number of genomes.
#' 
#' For every pair of genomes a result file is produced. If two genomes have genome_id's \samp{GID111},
#' and \samp{GID222} then the result file \samp{GID111_vs_GID222.txt} (or \samp{GID222_vs_GID111.txt}) will
#' be found in \samp{out.folder} after the completion of this search.
#' 
#' The \samp{out.folder} is scanned for already existing result files, and \code{\link{blastAllAll}} never
#' overwrites an existing result file. If a file with the name \samp{GID111_vs_GID222.txt} already exists in
#' the \samp{out.folder}, this particular search is skipped. This makes it possible to run multiple jobs in
#' parallell, writing to the same \samp{out.folder}. It also makes it possible to add new genomes, and only
#' BLAST the new combinations without repeating previous comparisons. 
#' 
#' This search can be slow if the genomes contain many proteins and it scales quadratically in the number of
#' input files. It is best suited for the study of a smaller number of genomes. By
#' starting multiple R sessions, you can speed up the search by running \code{\link{blastAllAll}} from each R
#' session, using the same \samp{out.folder} but different integers for the \code{job} option. If you are
#' using a computing cluster you can also increase the number of CPUs by increasing \code{threads}.
#' 
#' The result files are text files, and can be read into R using \code{\link{readBlastTable}}, but more
#' commonly they are used as input to \code{\link{bDist}} to compute distances between sequences for subsequent
#' clustering.
#' 
#' @return The function produces a result file for each pair of files listed in \samp{prot.files}N}.
#' These result files are located in \code{out.folder}. Existing files are never overwritten by
#' \code{\link{blastAllAll}}, if you want to re-compute something, delete the corresponding result files first.
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
#' @seealso \code{\link{panPrep}}, \code{\link{readBlastTable}}, \code{\link{bDist}}.
#' 
#' @examples 
#' \dontrun{
#' # This example requires the external BLAST+ software
#' # Using protein files in this package
#' xpth <- file.path(path.package("micropan"),"extdata")
#' pf <- file.path(xpth,c("Example_proteins_GID1.fasta.xz", "Example_proteins_GID2.fasta.xz", "Example_proteins_GID3.fasta.xz"))
#' 
#' # We need to uncompress them first...
#' prot.files <- tempfile(fileext = c("_GID1.fasta.xz","_GID2.fasta.xz","_GID3.fasta.xz"))
#' ok <- file.copy(from = pf, to = prot.files)
#' prot.files <- unlist(lapply(prot.files, xzuncompress))
#' 
#' # Blasting all versus all
#' out.dir <- "."
#' blastpAllAll(prot.files, out.folder = out.dir)
#' 
#' # Reading results, and computing blast.distances
#' blast.files <- list.files(out.dir, pattern="GID[0-9]+_vs_GID[0-9]+.txt")
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
blastpAllAll <- function(prot.files, out.folder, e.value = 1, job = 1, threads = 1, verbose = TRUE){
  if(available.external("blast+")){
    genome_id <- str_extract(prot.files, "GID[0-9]+")
    prot.files <- prot.files[order(genome_id)]
    genome_id <- genome_id[order(genome_id)]
    if(.Platform$OS.type == "windows"){
      outfmt <- "\"6 qseqid sseqid evalue bitscore\""
    } else {
      outfmt <- "'6 qseqid sseqid evalue bitscore'"
    }
    for(i in 1:length(prot.files)){
      log.fil <- file.path(out.folder, "log.txt")
      db.fil <- file.path(out.folder, str_c("blastDB", job))
      system(paste("makeblastdb -logfile",log.fil, "-dbtype prot -out", db.fil, "-in", prot.files[i]))
      for(j in i:length(prot.files)){
        rname <- str_c(genome_id[j], "_vs_", genome_id[i], ".txt")
        res.files <- list.files(out.folder, pattern = "txt$")
        if(!(rname %in% res.files)){
          if(verbose) cat("blastAllAll: ", rname, "\n")
          input <- str_c("-query ", prot.files[j])
          dbase <- str_c("-db ", db.fil)
          output <- str_c("-out ", file.path(out.folder, rname))
          cmd <- paste("blastp",
                       "-matrix BLOSUM45",
                       "-evalue", e.value,
                       "-num_threads", threads,
                       "-outfmt", outfmt,
                       input, dbase, output)
          system(cmd)
        }
      }
    }
    ok <- file.remove(str_c(db.fil, ".pin"))
    ok <- file.remove(str_c(db.fil, ".phr"))
    ok <- file.remove(str_c(db.fil, ".psq"))
    ok <- file.remove(log.fil, str_c(log.fil, ".perf"))
    invisible(TRUE)
  }
}
