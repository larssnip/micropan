#' @name panMatrix
#' @title Computing the pan-matrix for a set of gene clusters
#' 
#' @description A pan-matrix has one row for each genome and one column for each gene cluster, and
#' cell \samp{[i,j]} indicates how many members genome \samp{i} has in gene family \samp{j}.
#' 
#' @param clustering A vector of integers indicating the gene cluster for every sequence. Sequences
#' with the same number belong to the same cluster. The name of each element is the tag identifying
#' the sequence.
#' 
#' @details The pan-matrix is a central data structure for pan-genomic analysis. It is a matrix with
#' one row for each genome in the study, and one column for each gene cluster. Cell \samp{[i,j]}
#' contains an integer indicating how many members genome \samp{i} has in cluster \samp{j}.
#' 
#' The input \code{clustering} must be an integer vector with one element for each sequence in the study,
#' typically produced by either \code{\link{bClust}} or \code{\link{dClust}}. The name of each element
#' is a text identifying every sequence. The value of each element indicates the cluster, i.e. those
#' sequences with identical values are in the same cluster. IMPORTANT: The name of each sequence must
#' contain the GID-tag for each genome, i.e. they must of the form \samp{GID111_seq1}, \samp{GID111_seq2},...
#' where the \samp{GIDxxx} part indicates which genome the sequence belongs to. See \code{\link{panPrep}}
#' for details.
#' 
#' The rows of the pan-matrix is named by the GID-tag for every genome. The columns are just named
#' \samp{Cluster_x} where \samp{x} is an integer copied from \samp{clustering}.
#' 
#' @return An integer matrix with a row for each genome and a column for each sequence cluster.
#' The input vector \samp{clustering} is attached as the attribute \samp{clustering}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{bClust}}, \code{\link{dClust}}, \code{\link{distManhattan}},
#' \code{\link{distJaccard}}, \code{\link{fluidity}}, \code{\link{chao}},
#' \code{\link{binomixEstimate}}, \code{\link{heaps}}, \code{\link{rarefaction}}.
#' 
#' @examples 
#' # Loading clustering data in this package
#' data(xmpl.bclst)
#' 
#' # Pan-matrix based on the clustering
#' panmat <- panMatrix(xmpl.bclst)
#' 
#' # Plotting cluster distribution
#' tibble(Clusters = as.integer(table(factor(colSums(panmat > 0), levels = 1:nrow(panmat)))),
#'        Genomes = 1:nrow(panmat)) %>% 
#' ggplot(aes(x = Genomes, y = Clusters)) +
#' geom_col()
#' 
#' @importFrom stringr str_extract str_c
#' 
#' @export panMatrix
#' 
panMatrix <- function(clustering){
  gids <- str_extract(names(clustering), "GID[0-9]+")
  ugids <- sort(unique(gids))
  uclst <- sort(unique(clustering))
  pan.matrix <- matrix(0, nrow = length(ugids), ncol = length(uclst))
  rownames(pan.matrix) <- ugids
  colnames(pan.matrix) <- str_c("Cluster", uclst)
  for(i in 1:length(ugids)){
    tb <- table(clustering[gids == ugids[i]])
    idd <- as.numeric(names(tb))
    idx <- which(uclst %in% idd)
    pan.matrix[i,idx] <- tb
  }
  attr(pan.matrix, "clustering") <- clustering
  return(pan.matrix)
}


#' #' @rdname generic.Panmat
#' #' @name plot.Panmat
#' #' @title Plot and summary of \code{Panmat} objects
#' #' 
#' #' @description Generic functions for plotting and printing the content of a \code{Panmat} object.
#' #' 
#' #' @param x A \code{Panmat} object, see below.
#' #' @param object A \code{Panmat} object, see below.
#' #' @param col The color, default is \samp{"black"}, of interior and borders of the bars in the barplot.
#' #' @param xlab The label of the X axis.
#' #' @param ylab The label of the Y axis.
#' #' @param \dots Optional (graphical) arguments.
#' #' 
#' #' @details A \code{Panmat} object contains a pan-matrix, which is the fundamental data structure
#' #' for pan-genome analyses. It is a small (S3) extension to a \code{matrix}. It has one row for each
#' #' genome in the study, and one column for each gene cluster. The number in cell \samp{[i,j]} is the
#' #' number of sequences in genome \samp{i} that belongs to cluster \samp{j}. A \code{Panmat} object is
#' #' typically created by the function \code{\link{panMatrix}}.
#' #' 
#' #' The \code{\link{plot.Panmat}} function will display the content of the \code{Panmat} object as a bar
#' #' chart showing the number of clusters found in 1,2,...,G genomes, where G is the total number of genomes
#' #' in the study (rows in \samp{Panmat}).
#' #' 
#' #' The \code{\link{summary.Panmat}} function will display a text giving the same information as
#' #' \code{\link{plot.Panmat}}.
#' #' 
#' #' @author Lars Snipen and Kristian Hovde Liland.
#' #' 
#' #' @seealso \code{\link{panMatrix}}.
#' #' 
#' #' @examples # See examples in the Help-file for panMatrix.
#' #' 
#' #' @export
#' plotPanmatrix <- function(x, col="black", xlab="Number of genomes", ylab="Number of clusters", ... ){
#'   # x is a Panmat
#'   x[which( x > 0, arr.ind=T )] <- 1
#'   levs <- 1:dim( x )[1]
#'   y <- table( factor( colSums( x ), levels=levs ) )
#'   barplot( y, col=col, border=col, names.arg=levs, xlab=xlab, ylab=ylab, ... )
#' }
#' #' @rdname generic.Panmat
#' #' @export
#' summary.Panmat <- function( object, ... ){
#'   # object is a Panmat
#'   object[which( object > 0, arr.ind=T )] <- 1
#'   levs <- 1:nrow( object )
#'   y <- table( factor( colSums( object ), levels=levs ) )
#'   for( i in 1:length( levs ) ){
#'     cat( y[i], "clusters found in", levs[i], "genomes\n" )
#'   }
#' }


