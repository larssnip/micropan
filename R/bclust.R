#' @name bClust
#' @title Clustering sequences based on pairwise distances
#' 
#' @description Sequences are clustered by hierarchical clustering based on a set of pariwise distances.
#' The distances must take values between 0.0 and 1.0, and all pairs \emph{not} listed are assumed to
#' have distance 1.0.
#' 
#' @param dist.tbl A \code{tibble} with pairwise distances. 
#' @param linkage A text indicating what type of clustering to perform, either \samp{complete} (default),
#' \samp{average} or \samp{single}.
#' @param threshold Specifies the maximum size of a cluster. Must be a distance, i.e. a number between
#' 0.0 and 1.0.
#' @param verbose Logical, turns on/off text output during computations.
#' 
#' @details  Computing clusters (gene families) is an essential step in many comparative studies.
#' \code{\link{bClust}} will assign sequences into gene families by a hierarchical clustering approach.
#' Since the number of sequences may be huge, a full all-against-all distance matrix will be impossible
#' to handle in memory. However, most sequence pairs will have an \sQuote{infinite} distance between them,
#' and only the pairs with a finite (smallish) distance need to be considered.
#'
#' This function takes as input the distances in \code{dist.tbl} where only the relevant distances are
#' listed. The columns \samp{Query} and \samp{Hit} contain tags identifying pairs of sequences. The column
#' \samp{Distance} contains the distances, always a number from 0.0 to 1.0. Typically, this is the output
#' from \code{\link{bDist}}. All pairs of sequences \emph{not} listed are assumed to have distance 1.0,
#' which is considered the \sQuote{infinite} distance.
#' All sequences must be listed at least once in ceither column \samp{Query} or \samp{Hit} of the \code{dist.tbl}.
#' This should pose no problem, since all sequences must have distance 0.0 to themselves, and should be listed
#' with this distance once (\samp{Query} and \samp{Hit} containing the same tag). 
#' 
#' The \samp{linkage} defines the type of clusters produced. The \samp{threshold} indicates the size of
#' the clusters. A \samp{single} linkage clustering means all members of a cluster have at least one other
#' member of the same cluster within distance \samp{threshold} from itself. An \samp{average} linkage means
#' all members of a cluster are within the distance \samp{threshold} from the center of the cluster. A
#' \samp{complete} linkage means all members of a cluster are no more than the distance \samp{threshold}
#' away from any other member of the same cluster. 
#' 
#' Typically, \samp{single} linkage produces big clusters where members may differ a lot, since they are
#' only required to be close to something, which is close to something,...,which is close to some other
#' member. On the other extreme, \samp{complete} linkage will produce small and tight clusters, since all
#' must be similar to all. The \samp{average} linkage is between, but closer to \samp{complete} linkage. If
#' you want the \samp{threshold} to specify directly the maximum distance tolerated between two members of
#' the same gene family, you must use \samp{complete} linkage. The \samp{single} linkage is the fastest
#' alternative to compute. Using the default setting of \samp{single} linkage and maximum \samp{threshold}
#' (1.0) will produce the largest and fewest clusters possible.
#' 
#' @return The function returns a vector of integers, indicating the cluster membership of every unique
#' sequence from the \samp{Query} or \samp{Hit} columns of the input \samp{dist.tbl}. The name
#' of each element indicates the sequence. The numerical values have no meaning as such, they are simply
#' categorical indicators of cluster membership.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{bDist}}, \code{\link{hclust}}, \code{\link{dClust}}, \code{\link{isOrtholog}}.
#' 
#' @examples
#' # Loading example BLAST distances
#' data(xmpl.bdist)
#' 
#' # Clustering with default settings
#' clst <- bClust(xmpl.bdist)
#' # Other settings, and verbose
#' clst <- bClust(xmpl.bdist, linkage = "average", threshold = 0.5, verbose = TRUE)
#' 
#' @importFrom igraph graph_from_edgelist components degree
#' @importFrom stats hclust as.dist cutree
#' @importFrom dplyr filter %>% 
#' @importFrom rlang .data
#' 
#' @export bClust
#' 
bClust <- function(dist.tbl, linkage = "complete", threshold = 0.75, verbose = TRUE){
  if(verbose) cat("bClust:\n")
  linknum <- grep(linkage, c("single", "average", "complete"))
  dist.tbl %>% 
    filter(.data$Distance < threshold) -> dist.tbl
  utag <- sort(unique(c(dist.tbl$Dbase, dist.tbl$Query))) # Important to sort here!
    
  if(verbose) cat("...constructing graph with", length(utag), "sequences (nodes) and", nrow(dist.tbl), "distances (edges)\n")
  M <- matrix(as.numeric(factor(c(dist.tbl$Dbase, dist.tbl$Query), levels = utag)), ncol = 2, byrow = F)
  g <- graph_from_edgelist(M, directed = F)
  cls <- components(g)
  if(verbose) cat("...found", cls$no, "single linkage clusters\n")
  tibble(cluster = cls$membership,
         tag = utag) -> cls.tbl
  
  if(linknum > 1){
    ucls <- sort(unique(cls.tbl$cluster))
    incomplete <- which(sapply(ucls, function(j){
      v <- which(cls$membership == j)
      degg <- degree(g, v)
      return(min(degg) < (length(degg) + 1))
    }))
    if(verbose) cat("...found", length(incomplete), "incomplete clusters\n")
    if(length(incomplete) > 0){
      cls.tbl %>% 
        filter(.data$cluster %in% incomplete) %>% 
        mutate(cluster = .data$cluster * 1000) -> inc.tbl
      cls.tbl %>% 
        filter(!(.data$cluster %in% incomplete)) -> cls.tbl
      ucls.c <- unique(inc.tbl$cluster)
      for(i in 1:length(ucls.c)){   # for each incomplete cluster
        inc.tbl %>% 
          filter(.data$cluster == ucls.c[i]) -> tbl
        D <- matrix(1, nrow = nrow(tbl), ncol = nrow(tbl))
        rownames(D) <- colnames(D) <- tbl$tag
        dist.tbl %>% 
          filter(.data$Dbase %in% tbl$tag | .data$Query %in% tbl$tag) %>% 
          mutate(Dbase = factor(.data$Dbase, levels = tbl$tag)) %>% 
          mutate(Query = factor(.data$Query, levels = tbl$tag)) -> d.tbl
        M <- matrix(c(as.integer(d.tbl$Dbase), as.integer(d.tbl$Query)), ncol = 2, byrow = F)
        D[M] <- d.tbl$Distance
        D[M[,c(2,1)]] <- d.tbl$Distance
        if(linknum == 2){
          clst <- hclust(as.dist(D), method = "average")
        } else {
          clst <- hclust(as.dist(D), method = "complete")
        }
        tbl %>% 
          mutate(cluster = .data$cluster + cutree(clst, h = threshold)) %>% 
          bind_rows(cls.tbl) -> cls.tbl
        if(verbose) cat(i, "/", length(ucls.c), "\r")
      }
    }
  }
  clustering <- as.integer(factor(cls.tbl$cluster))  # to get values 1,2,3,...
  names(clustering) <- cls.tbl$tag
  if(verbose) cat("...ended with", length(unique(clustering)),
                  "clusters, largest cluster has",
                  max(table(clustering)), "members\n")
  return(sort(clustering))
}


#' @name isOrtholog
#' @title Identifies orthologs in gene clusters
#' 
#' @description Finds the ortholog sequences in every cluster based on pairwise distances.
#' 
#' @param clustering A vector of integers indicating the cluster for every sequence. Sequences with
#' the same number belong to the same cluster. The name of each element is the tag identifying the sequence.
#' @param dist.tbl A \code{tibble} with pairwise distances. The columns \samp{Query} and
#' \samp{Hit} contain tags identifying pairs of sequences. The column \samp{Distance} contains the
#' distances, always a number from 0.0 to 1.0.
#' 
#' @details The input \code{clustering} is typically produced by \code{\link{bClust}}. The input
#' \code{dist.tbl} is typically produced by \code{\link{bDist}}.
#' 
#' The concept of orthologs is difficult for prokaryotes, and this function finds orthologs in a
#' simplistic way. For a given cluster, with members from many genomes, there is one ortholog from every
#' genome. In cases where a genome has two or more members in the same cluster, only one of these is an
#' ortholog, the rest are paralogs.
#' 
#' Consider all sequences from the same genome belonging to the same cluster. The ortholog is defined as
#' the one having the smallest sum of distances to all other members of the same cluster, i.e. the one
#' closest to the \sQuote{center} of the cluster.
#' 
#' Note that the status as ortholog or paralog depends greatly on how clusters are defined in the first
#' place. If you allow large and diverse (and few) clusters, many sequences will be paralogs. If you define
#' tight and homogenous (and many) clusters, almost all sequences will be orthologs.
#' 
#' @return A vector of logicals with the same number of elements as the input \samp{clustering}, indicating
#' if the corresponding sequence is an ortholog (\code{TRUE}) or not (\code{FALSE}). The name of each
#' element is copied from \samp{clustering}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{bDist}}, \code{\link{bClust}}.
#' 
#' @examples
#' \dontrun{
#' # Loading distance data and their clustering results
#' data(list = c("xmpl.bdist","xmpl.bclst"))
#' 
#' # Finding orthologs
#' is.ortholog <- isOrtholog(xmpl.bclst, xmpl.bdist)
#' # The orthologs are
#' which(is.ortholog)
#' }
#' 
#' @export isOrtholog
#' 
isOrtholog <- function(clustering, dist.tbl){
  uclst <- unique(clustering)
  tags <- names(clustering)
  is.ortholog <- rep(F, length(clustering))
  names(is.ortholog) <- tags
  for(i in 1:length(uclst)){
    idx <- which(clustering == uclst[i])
    idd <- which((dist.tbl$Query %in% tags[idx]) & (dist.tbl$Hit %in% tags[idx]))
    gidz <- str_extract(tags[idx], "GID[0-9]+")
    if(max(table(gidz)) > 1){
      D <- matrix(1, nrow = length(idx), ncol = length(idx))
      a <- as.numeric(factor(dist.tbl$Query[idd], levels = tags[idx]))
      b <- as.numeric(factor(dist.tbl$Hit[idd], levels = tags[idx]))
      D[matrix(c(a,b), ncol = 2, byrow = F)] <- dist.tbl$Distance[idd]
      D[matrix(c(b,a), ncol = 2, byrow = F)] <- dist.tbl$Distance[idd]
      ixx <- order(rowSums(D))
      ixd <- which(!duplicated(gidz[ixx]))
      is.ortholog[idx[ixx[ixd]]] <- T
    } else {
      is.ortholog[idx] <- T
    }
    cat(i, "/", length(uclst), "\r")
  }
  return(is.ortholog)
}
   
    
    
