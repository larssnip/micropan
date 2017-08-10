#' @name mpneumoniae
#' @aliases Mpneumoniae Mpneumoniae.table Mpneumoniae.blast.distances Mpneumoniae.blast.clustering Mpneumoniae.blast.panmat Mpneumoniae.domain.clustering Mpneumoniae.domain.panmat
#' @docType data
#' @title Data sets for the \emph{Mycoplasma pneumoniae} casestudy
#' 
#' @description This data set contains several files with various objects related to the casestudy
#' example used for illustration purposes in the \code{micropan} package.
#' 
#' @usage 
#' data(Mpneumoniae.table)
#' data(Mpneumoniae.blast.distances)
#' data(Mpneumoniae.blast.clustering)
#' data(Mpneumoniae.blast.panmat)
#' data(Mpneumoniae.domain.clustering)
#' data(Mpneumoniae.domain.panmat)
#' 
#' @details 
#' \samp{Mpneumoniae.table} is a \code{data.frame} with 7 rows holding some information about
#' the 7 genomes in the casestudy.
#' 
#' \samp{Mpneumoniae.blast.distances} is a \code{data.frame} with 3 columns holding all
#' computed BLAST distances between pairs of sequences in the 7 genomes. This \code{data.frame}
#' has 139 543 rows.
#' 
#' \samp{Mpneumoniae.blast.clustering} is a clustering vector of all the 9573 sequences in the
#' genomes based on \samp{Mpneumoniae.blast.distances}.
#' 
#' \samp{Mpneumoniae.blast.panmat} is a \code{Panmat} object containing a pan-matrix with 7
#' rows and 1210 columns based on \samp{Mpneumoniae.blast.clustering}.
#' 
#' \samp{Mpneumoniae.domain.clustering} is a clustering vector of 5265 sequences in the genomes
#' based on domain sequences. Notice that only sequences having at least one protein domain is
#' considered here (5265 out of the total 9573).
#' 
#' \samp{Mpneumoniae.domain.panmat} is a \code{Panmat} object containing a pan-matrix with 7 rows
#' and 445 columns based on \samp{Mpneumoniae.domain.clustering}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples 
#' # Genome overview table
#' data(Mpneumoniae.table) #loads the Mpneumoniae.table
#' if(interactive()){
#'    View(Mpneumoniae.table)
#' } else {
#'    str(Mpneumoniae.table)
#' }
#' 
#' # BLAST distances, only the first 20 are displayed
#' data(Mpneumoniae.blast.distances) #loads the Mpneumoniae.blast.distances
#' if(interactive()){
#'    View(Mpneumoniae.blast.distances[1:20,])
#' } else {
#'    str(Mpneumoniae.blast.distances[1:20,])
#' }
#' 
#' # BLAST clustering vector
#' data(Mpneumoniae.blast.clustering) #loads the Mpneumoniae.blast.clustering
#' Mpneumoniae.blast.clustering[1:30]
#' 
#' # BLAST pan-matrix
#' data(Mpneumoniae.blast.panmat) #loads the Mpneumoniae.blast.panmat
#' summary(Mpneumoniae.blast.panmat)
#' 
#' # Domain sequence clustering vector
#' data(Mpneumoniae.domain.clustering) #loads the Mpneumoniae.domain.clustering
#' Mpneumoniae.domain.clustering[1:30]
#' 
#' # Domain sequence pan-matrix
#' data(Mpneumoniae.domain.panmat) #loads the Mpneumoniae.domain.panmat
#' summary(Mpneumoniae.domain.panmat)
#' 
NULL