#' @name xmpl
#' @aliases xmpl xmpl.bdist xmpl.bclst xmpl.panmat
#' @docType data
#' @title Data sets for use in examples
#' 
#' @description This data set contains several files with various objects used in examples
#' in some of the functions in the \code{micropan} package.
#' 
#' @usage 
#' data(xmpl.bdist)
#' data(xmpl.bclst)
#' data(xmpl.panmat)
#' 
#' @details 
#' \samp{xmpl.bdist} is a \code{tibble} with 4 columns holding all
#' BLAST distances between pairs of proteins in an example with 10 small genomes.
#' 
#' \samp{xmpl.bclst} is a clustering vector of all proteins in the
#' genomes from \samp{xmpl.bdist}.
#' 
#' \samp{xmpl.panmat} is a pan-matrix with 10 rows and 1210 columns
#' computed from \samp{xmpl.bclst}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @examples 
#' 
#' # BLAST distances, only the first 20 are displayed
#' data(xmpl.bdist)
#' head(xmpl.bdist)
#' 
#' # Clustering vector
#' data(xmpl.bclst)
#' print(xmpl.bclst[1:30])
#' 
#' # Pan-matrix
#' data(xmpl.panmat)
#' head(xmpl.panmat)
#' 
NULL