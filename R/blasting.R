#' @name blastAllAll
#' @title Making BLAST search all against all genomes
#' 
#' @description Runs a reciprocal all-against-all BLAST search to look for similarity of proteins
#' within and across genomes. The main job is done by the BLAST+ software.
#' 
#' @param prot.files A text vector with the names of the FASTA files where the protein sequences of
#' each genome is found.
#' @param out.folder The name of the folder where the result files should end up.
#' @param e.value The chosen E-value threshold in BLAST. Default is \samp{e.value=1}, a smaller value
#' will speed up the search at the cost of less sensitivity.
#' @param job An integer to separate multiple jobs. You may want to run several jobs in parallell,
#' and each job should have different number here to avoid confusion on databases. Default is
#' \samp{job=1}.
#' @param threads The number of CPU's to use.
#' @param verbose Logical, if \code{TRUE} some text output is produced to monitor the progress.
#' 
#' @details A basic step in pangenomics and many other comparative studies is to cluster proteins into
#' groups or families. One commonly used approach is based on reciprocal BLASTing. This function uses the
#' BLAST+ software available for free from NCBI (Camacho et al, 2009). 
#' 
#' A vector listing FASTA files of protein sequences is given as input in \samp{prot.files}. These files
#' must have the GID-tag in the first token of every header, and in their filenames as well, i.e. all input
#' files should first be prepared by \code{\link{panPrep}} to ensure this. Note that only protein sequences
#' are considered here. If your coding genes are stored as DNA, please translate them to protein prior to
#' using this function, see \code{\link[microseq]{translate}}.
#' 
#' A BLAST database is made from each genome in turn. Then all genomes are queried against this database,
#' and for every pair of genomes a result file is produced. If two genomes have GID-tags \samp{GID111},
#' and \samp{GID222} then both result file \samp{GID111_vs_GID222.txt} and \samp{GID222_vs_GID111.txt} will
#' be found in \samp{out.folder} after the completion of this search. This reciprocal (two-way) search is
#' required because of the heuristics of BLAST.
#' 
#' The \samp{out.folder} is scanned for already existing result files, and \code{\link{blastAllAll}} never
#' overwrites an existing result file. If a file with the name \samp{GID111_vs_GID222.txt} already exists in
#' the \samp{out.folder}, this particular search is skipped. This makes it possible to run multiple jobs in
#' parallell, writing to the same \samp{out.folder}. It also makes it possible to add new genomes, and only
#' BLAST the new combinations without repeating previous comparisons. 
#' 
#' This search can be slow if the genomes contain many proteins and it scales quadratically in the number of
#' input files. It is best suited for the study of a smaller number of genomes (less than say 100). By
#' starting multiple R sessions, you can speed up the search by running \code{\link{blastAllAll}} from each R
#' session, using the same \samp{out.folder} but different integers for the \code{job} option. If you are
#' using a computing cluster you can also increase the number of CPUs by increasing \code{threads}.
#' 
#' The result files are text files, and can be read into R using \code{\link{readBlastTable}}, but more
#' commonly they are used as input to \code{\link{bDist}} to compute distances between sequences for subsequent
#' clustering.
#' 
#' @return The function produces \emph{N*N} result files if \samp{prot.files} lists \emph{N} sequence files.
#' These result files are located in \code{out.folder}. Existing files are never overwritten by
#' \code{\link{blastAllAll}}, if you want to re-compute something, delete the corresponding result files first.
#' 
#' @references Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L.
#' (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10:421.
#' 
#' @note The BLAST+ software must be installed on the system for this function to work, i.e. the command
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
#' prot.files <- file.path(xpth,c("Example_proteins_GID1.fasta.xz",
#'                                "Example_proteins_GID2.fasta.xz",
#'                                "Example_proteins_GID3.fasta.xz"))
#' 
#' # We need to uncompress them first...
#' tf <- tempfile(fileext=c("GID1.fasta.xz","GID2.fasta.xz","GID3.fasta.xz"))
#' s <- file.copy(prot.files,tf)
#' tf <- unlist(lapply(tf,xzuncompress))
#' 
#' # Blasting all versus all...(requires BLAST+)
#' tmp.dir <- tempdir()
#' blastAllAll(tf,out.folder=tmp.dir)
#' 
#' # Reading results, and computing blast.distances
#' blast.files <- dir(tmp.dir,pattern="GID[0-9]+_vs_GID[0-9]+.txt")
#' blast.distances <- bDist(file.path(tmp.dir,blast.files))
#' 
#' # ...and cleaning tmp.dir...
#' s <- file.remove(tf)
#' s <- file.remove(file.path(tmp.dir,blast.files))
#' }
#' 
#' @importFrom microseq gregexpr
#' @export
blastAllAll <- function(prot.files, out.folder, e.value = 1, job = 1, threads = 1, verbose = T){
  if(available.external("blast+")){
    for(i in 1:length(prot.files)){
      log.fil <- file.path(out.folder, "log.txt")
      db.fil <- file.path(out.folder, paste0("blastDB", job))
      command <- paste("makeblastdb -logfile",log.fil, "-dbtype prot -out", db.fil, "-in", prot.files[i])
      system(command)
      gi <- gregexpr("GID[0-9]+", prot.files[i], extract = T)
      for(j in 1:length(prot.files)){
        gj <- gregexpr("GID[0-9]+", prot.files[j], extract = T)    
        rname <- paste0(gj, "_vs_", gi, ".txt")
        res.files <- dir(out.folder)
        if(!(rname %in% res.files)){
          if(verbose) cat("blastAllAll: ", rname, "\n")
          input <- paste0("-query ", prot.files[j])
          dbase <- paste0("-db ", db.fil)
          output <- paste0("-out ", file.path(out.folder, rname))
          command <- paste("blastp -matrix BLOSUM45 -evalue", e.value, "-num_threads", threads,
                            "-outfmt 6", input, dbase, output)
          system(command)
        }
      }
    }
    file.remove(paste0(db.fil, ".pin"))
    file.remove(paste0(db.fil, ".phr"))
    file.remove(paste0(db.fil, ".psq"))
    file.remove(log.fil, paste0(log.fil, ".perf"))
    return(TRUE)
  }
}


#' @name readBlastTable
#' @title Reading BLAST result file
#' 
#' @description Reading a file produced by the BLAST+ software set up to produce tabular output.
#' 
#' @param blast.file Name of file to read.
#' 
#' @details This function will read files produced by the BLAST+ software where the option \samp{-outfmt 6}
#' has been invoked during its call. This option forces BLAST to produce a short tabular text output for
#' each BLAST search. The function \code{\link{blastAllAll}} produces such files.
#' 
#' @return The content of the file is returned as a \samp{data.frame} with 12 columns and one row for each
#' BLAST result. The columns have self-explanatory names.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{blastAllAll}}, \code{\link{bDist}}.
#' 
#' @examples 
#' # Using a BLAST result file in this package
#' xpth <- file.path(path.package("micropan"),"extdata")
#' blast.file <- file.path( xpth, "GID1_vs_GID2.txt.xz" )
#' 
#' # We need to uncompress it first...
#' tf <- tempfile(fileext=".xz")
#' s <- file.copy(blast.file,tf)
#' tf <- xzuncompress(tf)
#' 
#' #...then we can read it...
#' blast.table <- readBlastTable(tf)
#' 
#' # ...and deleting temporary file
#' s <- file.remove(tf)
#' 
#' @importFrom utils read.table
#' 
#' @export
readBlastTable <- function(blast.file){
  columnames <- c("Query", "Hit", "Percent.identity", "Alignment.length", "Mismatches",
                  "Gap.openings", "Query.start", "Query.end", "Hit.start", "Hit.end", "E.value", "Bit.score")
  data <- read.table(blast.file, sep = "\t", quote = "", header = F, col.names = columnames, stringsAsFactors = F)
  return(data)
}






