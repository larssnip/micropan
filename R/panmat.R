#' @name panMatrix
#' @title Computing the pan-matrix for a set of gene clusters
#' 
#' @description A pan-matrix has one row for each genome and one column for each gene cluster, and
#' cell \samp{[i,j]} indicates how many members genome \samp{i} has in gene family \samp{j}.
#' 
#' @param clustering A named vector of integers.
#' 
#' @details The pan-matrix is a central data structure for pan-genomic analysis. It is a matrix with
#' one row for each genome in the study, and one column for each gene cluster. Cell \samp{[i,j]}
#' contains an integer indicating how many members genome \samp{i} has in cluster \samp{j}.
#' 
#' The input \code{clustering} must be a named integer vector with one element for each sequence in the study,
#' typically produced by either \code{\link{bClust}} or \code{\link{dClust}}. The name of each element
#' is a text identifying every sequence. The value of each element indicates the cluster, i.e. those
#' sequences with identical values are in the same cluster. IMPORTANT: The name of each sequence must
#' contain the \samp{genome_id} for each genome, i.e. they must of the form \samp{GID111_seq1}, \samp{GID111_seq2},...
#' where the \samp{GIDxxx} part indicates which genome the sequence belongs to. See \code{\link{panPrep}}
#' for details.
#' 
#' The rows of the pan-matrix is named by the \samp{genome_id} for every genome. The columns are just named
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
#' \dontrun{
#' # Plotting cluster distribution
#' library(ggplot2)
#' tibble(Clusters = as.integer(table(factor(colSums(panmat > 0), levels = 1:nrow(panmat)))),
#'        Genomes = 1:nrow(panmat)) %>% 
#' ggplot(aes(x = Genomes, y = Clusters)) +
#' geom_col()
#' }
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

