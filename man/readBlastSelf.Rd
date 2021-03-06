% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdist.R
\name{readBlastSelf}
\alias{readBlastSelf}
\alias{readBlastPair}
\title{Reads BLAST result files}
\usage{
readBlastSelf(blast.files, e.value = 1, verbose = TRUE)
}
\arguments{
\item{blast.files}{A text vector of filenames.}

\item{e.value}{A threshold E-value to immediately discard (very) poor BLAST alignments.}

\item{verbose}{Logical, indicating if textual output should be given to monitor the progress.}
}
\value{
The functions returns a table with columns \samp{Dbase}, \samp{Query}, \samp{Bitscore}
and \samp{Distance}. Each row corresponds to a pair of sequences (a Dbase and a Query sequence) having at least
one BLAST hit between
them. All pairs \emph{not} listed have distance 1.0 between them. You should normally bind the output from 
\code{readBlastSelf} to the ouptut from \code{readBlastPair} and use the result as input to \code{\link{bDist}}.
}
\description{
Reads files from a search with blastpAllAll
}
\details{
The filenames given as input must refer to BLAST result files produced by \code{\link{blastpAllAll}}.

With \code{readBlastSelf} you only read the self-alignment results, i.e. blasting a genome against itself. With
\code{readBlastPair} you read all the other files, i.e. different genomes compared. You may use all blast file
names as input to both, they will select the proper files based on their names, e.g. GID1_vs_GID1.txt is read
by \code{readBlastSelf} while GID2_vs_GID1.txt is read by \code{readBlastPair}.

Setting a small \samp{e.value} threshold will filter the alignment, and may speed up this and later processing,
but you may also loose some important alignments for short sequences.

Both these functions are used by \code{\link{bDist}}. The reason we provide them separately is to allow the user
to complete this file reading before calling \code{\link{bDist}}. If you have a huge number of files, a
skilled user may utilize parallell processing to speed up the reading. For normal size data sets (e.g. less than 100 genomes)
you should probably use \code{\link{bDist}} directly.
}
\examples{
# Using BLAST result files in this package...
prefix <- c("GID1_vs_GID1_",
            "GID2_vs_GID1_",
            "GID3_vs_GID1_",
            "GID2_vs_GID2_",
            "GID3_vs_GID2_",
            "GID3_vs_GID3_")
bf <- file.path(path.package("micropan"), "extdata", str_c(prefix, ".txt.xz"))

# We need to uncompress them first...
blast.files <- tempfile(pattern = prefix, fileext = ".txt.xz")
ok <- file.copy(from = bf, to = blast.files)
blast.files <- unlist(lapply(blast.files, xzuncompress))

# Reading self-alignment files, then the other files
self.tbl <- readBlastSelf(blast.files)
pair.tbl <- readBlastPair(blast.files)

# ...and cleaning...
ok <- file.remove(blast.files)

# See also examples for bDist

}
\seealso{
\code{\link{bDist}}, \code{\link{blastpAllAll}}.
}
\author{
Lars Snipen.
}
