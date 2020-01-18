#' @name rarefaction
#' @title Rarefaction curves for a pan-genome
#' 
#' @description Computes rarefaction curves for a number of random permutations of genomes.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param n.perm The number of random genome orderings to use. If \samp{n.perm=1} the fixed order of
#' the genomes in \samp{pan.matrix} is used.
#' 
#' @details A rarefaction curve is simply the cumulative number of unique gene clusters we observe as
#' more and more genomes are being considered. The shape of this curve will depend on the order of the
#' genomes. This function will typically compute rarefaction curves for a number of (\samp{n.perm})
#' orderings. By using a large number of permutations, and then averaging over the results, the effect
#' of any particular ordering is smoothed.
#' 
#' The averaged curve illustrates how many new gene clusters we observe for each new genome. If this
#' levels out and becomes flat, it means we expect few, if any, new gene clusters by sequencing more
#' genomes. The function \code{\link{heaps}} can be used to estimate population openness based on this
#' principle.
#' 
#' @return A table with the curves in the columns. The first column is the number of genomes, while 
#' all other columns are the cumulative number of clusters, one column for each permutation.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{heaps}}, \code{\link{panMatrix}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package 
#' data(xmpl.panmat)
#' 
#' # Rarefaction
#' rar.tbl <- rarefaction(xmpl.panmat, n.perm = 1000)
#' 
#' \dontrun{
#' # Plotting
#' library(ggplot2)
#' library(tidyr)
#' rar.tbl %>% 
#'   gather(key = "Permutation", value = "Clusters", -Genome) %>% 
#'   ggplot(aes(x = Genome, y = Clusters, group = Permutation)) +
#'     geom_line()
#' }
#' 
#' @importFrom dplyr %>% bind_cols
#' @importFrom tibble tibble as_tibble
#' 
#' @export rarefaction
#' 
rarefaction <- function(pan.matrix, n.perm = 1){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  nmat <- matrix(0, nrow = nrow(pan.matrix), ncol = n.perm)
  cm <- apply(pan.matrix, 2, cumsum)
  nmat[,1] <- rowSums(cm > 0)
  if(n.perm > 1){
    for(i in 2:n.perm){
      cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)
      nmat[,i] <- rowSums(cm > 0)
      cat(i, "/", n.perm, "\r")
    }
  }
  nmat <- rbind(rep(0, ncol(nmat)), nmat)
  tibble(Genomes = 0:nrow(pan.matrix)) %>% 
    bind_cols(as_tibble(nmat, .name_repair = "minimal")) -> rtbl
  colnames(rtbl) <- c("Genome", str_c("Perm", 1:n.perm))
  return(rtbl)
}

