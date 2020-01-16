#' @name fluidity
#' @title Computing genomic fluidity for a pan-genome
#' 
#' @description Computes the genomic fluidity, which is a measure of population diversity.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param n.sim An integer specifying the number of random samples to use in the computations.
#' 
#' @details  The genomic fluidity between two genomes is defined as the number of unique gene
#' families divided by the total number of gene families (Kislyuk et al, 2011). This is averaged
#' over \samp{n.sim} random pairs of genomes to obtain a population estimate.
#' 
#' The genomic fluidity between two genomes describes their degree of overlap with respect to gene
#' cluster content. If the fluidity is 0.0, the two genomes contain identical gene clusters. If it
#' is 1.0 the two genomes are non-overlapping. The difference between a Jaccard distance (see
#' \code{\link{distJaccard}}) and genomic fluidity is small, they both measure overlap between
#' genomes, but fluidity is computed for the population by averaging over many pairs, while Jaccard
#' distances are computed for every pair. Note that only presence/absence of gene clusters are
#' considered, not multiple occurrences.
#' 
#' The input \samp{pan.matrix} is typically constructed by \code{\link{panMatrix}}.
#' 
#' @return A vector with two elements, the mean fluidity and its sample standard deviation over
#' the \samp{n.sim} computed values.
#' 
#' @references Kislyuk, A.O., Haegeman, B., Bergman, N.H., Weitz, J.S. (2011). Genomic fluidity:
#' an integrative view of gene diversity within microbial populations. BMC Genomics, 12:32.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distJaccard}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package
#' data(xmpl.panmat)
#' 
#' # Fluidity based on this pan-matrix
#' fluid <- fluidity(xmpl.panmat)
#' 
#' @importFrom stats sd
#' 
#' @export fluidity
#' 
fluidity <- function(pan.matrix, n.sim = 10){
  pan.matrix[which(pan.matrix > 0, arr.ind=T)] <- 1
  flu <- rep(0, n.sim)
  for(i in 1:n.sim){
    ii <- sample(nrow(pan.matrix), 2)
    flu[i] <- (sum(pan.matrix[ii[1],] > 0 & pan.matrix[ii[2],] == 0)
               + sum(pan.matrix[ii[1],] == 0 & pan.matrix[ii[2],] > 0)) / (sum(pan.matrix[ii[1],]) + sum(pan.matrix[ii[2],]))
  }
  flu.vec <- c(Mean = mean(flu), Std = sd(flu))
  return(flu.vec)
}


#' @name distJaccard
#' @title Computing Jaccard distances between genomes
#' 
#' @description Computes the Jaccard distances between all pairs of genomes.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' 
#' @details The Jaccard index between two sets is defined as the size of the intersection of
#' the sets divided by the size of the union. The Jaccard distance is simply 1 minus the Jaccard index.
#' 
#' The Jaccard distance between two genomes describes their degree of overlap with respect to gene
#' cluster content. If the Jaccard distance is 0.0, the two genomes contain identical gene clusters.
#' If it is 1.0 the two genomes are non-overlapping. The difference between a genomic fluidity (see
#' \code{\link{fluidity}}) and a Jaccard distance is small, they both measure overlap between genomes,
#' but fluidity is computed for the population by averaging over many pairs, while Jaccard distances are
#' computed for every pair. Note that only presence/absence of gene clusters are considered, not multiple
#' occurrences.
#' 
#' The input \samp{pan.matrix} is typically constructed by \code{\link{panMatrix}}.
#' 
#' @return A \code{dist} object (see \code{\link{dist}}) containing all pairwise Jaccard distances
#' between genomes.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{fluidity}}, \code{\link{dist}}.
#' 
#' @examples
#' # Loading a pan-matrix in this package
#' data(xmpl.panmat)
#' 
#' # Jaccard distances
#' Jdist <- distJaccard(xmpl.panmat)
#' 
#' # Making a dendrogram based on the distances,
#' # see example for distManhattan
#' 
#' @export distJaccard
#' 
distJaccard <- function(pan.matrix){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  D <- matrix(0, nrow = nrow(pan.matrix), ncol = nrow(pan.matrix))
  rownames(D) <- colnames(D) <- rownames(pan.matrix)
  for(i in 1:(nrow(pan.matrix) - 1)){
    for(j in (i+1):nrow(pan.matrix)){
      cs <- pan.matrix[i,] + pan.matrix[j,]
      D[j,i] <- D[i,j] <- 1 - sum(cs > 1)/sum(cs > 0)
    }
  }
  return(as.dist(D))
}


#' @name distManhattan
#' @title Computing Manhattan distances between genomes
#' 
#' @description Computes the (weighted) Manhattan distances beween all pairs of genomes.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param scale An optional scale to control how copy numbers should affect the distances.
#' @param weights Vector of optional weights of gene clusters.
#' 
#' @details The Manhattan distance is defined as the sum of absolute elementwise differences between
#' two vectors. Each genome is represented as a vector (row) of integers in \samp{pan.matrix}. The
#' Manhattan distance between two genomes is the sum of absolute difference between these rows. If
#' two rows (genomes) of the \samp{pan.matrix} are identical, the corresponding Manhattan distance
#' is \samp{0.0}.
#' 
#' The \samp{scale} can be used to control how copy number differences play a role in the distances
#' computed. Usually we assume that going from 0 to 1 copy of a gene is the big change of the genome,
#' and going from 1 to 2 (or more) copies is less. Prior to computing the Manhattan distance, the
#' \samp{pan.matrix} is transformed according to the following affine mapping: If the original value in
#' \samp{pan.matrix} is \samp{x}, and \samp{x} is not 0, then the transformed value is \samp{1 + (x-1)*scale}.
#' Note that with \samp{scale=0.0} (default) this will result in 1 regardless of how large \samp{x} was.
#' In this case the Manhattan distance only distinguish between presence and absence of gene clusters.
#' If \samp{scale=1.0} the value \samp{x} is left untransformed. In this case the difference between 1
#' copy and 2 copies is just as big as between 1 copy and 0 copies. For any \samp{scale} between 0.0 and
#' 1.0 the transformed value is shrunk towards 1, but a certain effect of larger copy numbers is still
#' present. In this way you can decide if the distances between genomes should be affected, and to what
#' degree, by differences in copy numbers beyond 1. Notice that as long as \samp{scale=0.0} (and no
#' weighting) the Manhattan distance has a nice interpretation, namely the number of gene clusters that
#' differ in present/absent status between two genomes.
#' 
#' When summing the difference across gene clusters we can also up- or downweight some clusters compared
#' to others. The vector \samp{weights} must contain one value for each column in \samp{pan.matrix}. The
#' default is to use flat weights, i.e. all clusters count equal. See \code{\link{geneWeights}} for
#' alternative weighting strategies.
#' 
#' @return A \code{dist} object (see \code{\link{dist}}) containing all pairwise Manhattan distances
#' between genomes.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distJaccard}}, \code{\link{geneWeights}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package
#' data(xmpl.panmat)
#' 
#' # Manhattan distances between genomes
#' Mdist <- distManhattan(xmpl.panmat)
#' 
#' \dontrun{
#' # Making a dendrogram based on shell-weighted distances
#' library(ggdendro)
#' weights <- geneWeights(xmpl.panmat, type = "shell")
#' Mdist <- distManhattan(xmpl.panmat, weights = weights)
#' ggdendrogram(dendro_data(hclust(Mdist, method = "average")),
#'   rotate = TRUE, theme_dendro = FALSE) +
#'   labs(x = "Genomes", y = "Shell-weighted Manhattan distance", title = "Pan-genome dendrogram")
#' }
#' 
#' @importFrom stats dist
#' 
#' @export distManhattan
#' 
distManhattan <- function(pan.matrix, scale = 0.0, weights = rep(1, ncol(pan.matrix))){
  if((scale > 1) | (scale < 0)){
    warning( "scale should be between 0.0 and 1.0, using scale = 0.0" )
    scale <- 0.0
  }
  idx <- which(pan.matrix > 0, arr.ind = T)
  pan.matrix[idx] <- 1 + (pan.matrix[idx] - 1) * scale
  pan.matrix <- t(t(pan.matrix) * weights)
  return(dist(pan.matrix, method = "manhattan"))
}


#' @name geneWeights
#' @title Gene cluster weighting
#' 
#' @description This function computes weights for gene cluster according to their distribution in a pan-genome.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param type A text indicating the weighting strategy.
#' 
#' @details When computing distances between genomes or a PCA, it is possible to give weights to the
#' different gene clusters, emphasizing certain aspects.
#' 
#' As proposed by Snipen & Ussery (2010), we have implemented two types of weighting: The default
#' \samp{"shell"} type means gene families occuring frequently in the genomes, denoted shell-genes, are
#' given large weight (close to 1) while those occurring rarely are given small weight (close to 0).
#' The opposite is the \samp{"cloud"} type of weighting. Genes observed in a minority of the genomes are
#' referred to as cloud-genes. Presumeably, the \samp{"shell"} weighting will give distances/PCA reflecting
#' a more long-term evolution, since emphasis is put on genes who have just barely diverged away from the
#' core. The \samp{"cloud"} weighting emphasizes those gene clusters seen rarely. Genomes with similar
#' patterns among these genes may have common recent history. A \samp{"cloud"} weighting typically gives
#' a more erratic or \sQuote{noisy} picture than the \samp{"shell"} weighting.
#' 
#' @return A vector of weights, one for each column in \code{pan.matrix}.
#' 
#' @references Snipen, L., Ussery, D.W. (2010). Standard operating procedure for computing pangenome
#' trees. Standards in Genomic Sciences, 2:135-141.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distManhattan}}.
#' 
#' @examples 
#' # See examples for distManhattan
#' 
#' @export geneWeights
#' 
geneWeights <- function(pan.matrix, type = c("shell", "cloud")){
  ng <- dim( pan.matrix )[1]
  nf <- dim( pan.matrix )[2]
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  cs <- colSums(pan.matrix)
  
  midx <- grep(type[1], c("shell", "cloud"))
  if(length(midx) == 0){
    warning("Unknown weighting:", type, ", using shell weights")
    midx <- 1
  }
  W <- rep(1, ncol(pan.matrix))
  x <- 1:nrow(pan.matrix)
  ww <- 1 / (1 + exp(((x - 1) - (max(x) - 1)/2) / ((max(x) - 1) / 10)))
  if(midx == 1) ww <- 1 - ww
  for(i in x) W[cs == i] <- ww[i]
  return(W)
}
