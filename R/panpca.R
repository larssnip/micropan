#' @name panPca
#' @title Principal component analysis of a pan-matrix
#' 
#' @description Computes a principal component decomposition of a pan-matrix, with possible
#' scaling and weightings.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param scale An optional scale to control how copy numbers should affect the distances.
#' @param weights Vector of optional weights of gene clusters.
#' 
#' @details A principal component analysis (PCA) can be computed for any matrix, also a pan-matrix.
#' The principal components will in this case be linear combinations of the gene clusters. One major
#' idea behind PCA is to truncate the space, e.g. instead of considering the genomes as points in a
#' high-dimensional space spanned by all gene clusters, we look for a few \sQuote{smart} combinations
#' of the gene clusters, and visualize the genomes in a low-dimensional space spanned by these directions.
#' 
#' The \samp{scale} can be used to control how copy number differences play a role in the PCA. Usually
#' we assume that going from 0 to 1 copy of a gene is the big change of the genome, and going from 1 to
#' 2 (or more) copies is less. Prior to computing the PCA, the \samp{pan.matrix} is transformed according
#' to the following affine mapping: If the original value in \samp{pan.matrix} is \samp{x}, and \samp{x}
#' is not 0, then the transformed value is \samp{1 + (x-1)*scale}. Note that with \samp{scale=0.0}
#' (default) this will result in 1 regardless of how large \samp{x} was. In this case the PCA only
#' distinguish between presence and absence of gene clusters. If \samp{scale=1.0} the value \samp{x} is
#' left untransformed. In this case the difference between 1 copy and 2 copies is just as big as between
#' 1 copy and 0 copies. For any \samp{scale} between 0.0 and 1.0 the transformed value is shrunk towards
#' 1, but a certain effect of larger copy numbers is still present. In this way you can decide if the PCA
#' should be affected, and to what degree, by differences in copy numbers beyond 1.
#' 
#' The PCA may also up- or downweight some clusters compared to others. The vector \samp{weights} must
#' contain one value for each column in \samp{pan.matrix}. The default is to use flat weights, i.e. all
#' clusters count equal. See \code{\link{geneWeights}} for alternative weighting strategies.
#' 
#' @return A \code{list} with three tables:
#' 
#' \samp{Evar.tbl} has two columns, one listing the component number and one listing the relative 
#' explained variance for each component. The relative explained variance always sums to 1.0 over
#' all components. This value indicates the importance of each component, and it is always in
#' descending order, the first component being the most important.
#' This is typically the first result you look at after a PCA has been computed, as it indicates
#' how many components (directions) you need to capture the bulk of the total variation in the data.
#' 
#' \samp{Scores.tbl} has a column listing the \samp{GID.tag} for each genome, and then one column for each
#' principal component. The columns are ordered corresponding to the elements in \samp{Evar}. The
#' scores are the coordinates of each genome in the principal component space.
#' 
#' \samp{Loadings.tbl} is similar to \samp{Scores.tbl} but contain values for each gene cluster
#' instead of each genome. The columns are ordered corresponding to the elements in \samp{Evar}.
#' The loadings are the contributions from each gene cluster to the principal component directions.
#' NOTE: Only gene clusters having a non-zero variance is used in a PCA. Gene clusters with the
#' same value for every genome have no impact and are discarded from the \samp{Loadings}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{distManhattan}}, \code{\link{geneWeights}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package
#' data(xmpl.panmat)
#' 
#' # Computing panPca
#' ppca <- panPca(xmpl.panmat)
#' 
#' # Plotting explained variance
#' ggplot(ppca$Evar.tbl) +
#'   geom_col(aes(x = Component, y = Explained.variance))
#' # Plotting scores
#' ggplot(ppca$Scores.tbl) +
#'   geom_text(aes(x = PC1, y = PC2, label = GID.tag))
#' # Plotting loadings
#' ggplot(ppca$Loadings.tbl) +
#'   geom_text(aes(x = PC1, y = PC2, label = Cluster))
#' 
#' @importFrom tibble as_tibble tibble
#' 
#' @export panPca
#' 
panPca <- function(pan.matrix, scale = 0.0, weights = rep(1, ncol(pan.matrix))){
  if((scale > 1) | (scale < 0)){
    warning("scale should be between 0.0 and 1.0, using scale=0.0")
    scale <- 0.0
  }
  idx <- which(pan.matrix > 0, arr.ind = T)
  pan.matrix[idx] <- 1 + (pan.matrix[idx] - 1) * scale
  pan.matrix <- t(t(pan.matrix) * weights)
  X <- pan.matrix[,which(apply(pan.matrix, 2, sd) > 0)]
  pca <- prcomp(X)
  pca.lst <- list(Evar.tbl     = tibble(Component = 1:length(pca$sdev),
                                        Explained.variance = pca$sdev^2/sum(pca$sdev^2)),
                  Scores.tbl   = as_tibble(pca$x, rownames = "GID.tag"),
                  Loadings.tbl = as_tibble(pca$rotation, rownames = "Cluster"))
  return(pca.lst)
}



#' #' @rdname scores.loadings
#' #' @name plotScores
#' #' @title Plotting scores and loadings in a \code{pca} object
#' #' 
#' #' @description Creates informative plots for a principal component analysis of a pan-matrix.
#' #' 
#' #' @param pca A \code{pca} object, see \code{\link{panPca}} for details.
#' #' @param x The component to display along the horizontal axis.
#' #' @param y The component to display along the vertical axis.
#' #' @param show.labels Logical indicating if labels should be displayed.
#' #' @param labels Alternative labels to use in the score-plot, see below.
#' #' @param col Colors for the points/labels, see below.
#' #' @param pch Marker type, see \code{\link{points}}.
#' #' @param \dots Additional arguments passed on to \code{points} or \code{text} (if labels are specified).
#' #' 
#' #' @details A \code{Panpca} object contains the results of a principal component analysis on a pan-matrix,
#' #' see \code{\link{panpca}} for details.
#' #' 
#' #' The \code{\link{plotScores}} gives a visual overview of how the genomes are positioned relative to
#' #' each other in the pan-genome space. The score-matrix of a \code{Panpca} has one row for each genome.
#' #' The original pan-matrix also has one row for each genome. Two genomes can be compared by their
#' #' corresponding rows in the pan-matrix, but can also be compared by their rows in the score-matrix,
#' #' and the latter matrix has (much) fewer columns designed to contain maximum of the original data
#' #' variation. A plot of the scores will give an approximate overview of how the genomes are located
#' #' relative to each other.
#' #' 
#' #' The \code{\link{plotLoadings}} gives a visual overview of how the gene clusters affect the principal
#' #' components. The loadings is a matrix with one row for each of the original non-core gene clusters
#' #' (core gene clusters have no variation across genomes). Clusters located close to the origin have
#' #' little impact. Clusters far from the origin has high impact, indicating they separate groups of genomes.
#' #' 
#' #' These two plots together can reveal information about the pan-genome: The score-plot shows if genomes
#' #' are grouped/separated, and the loading-plot can then tell you which gene clusters have high impact on
#' #' this grouping/separation.
#' #' 
#' #' The arguments \samp{x} and \samp{y} can be used to plot other components than component 1 and 2
#' #' (which is always the most informative). In some cases more components are needed to establish a
#' #' good picture, i.e. the explained variance is low for component 1 and 2 (see \code{\link{plot.Panpca}}
#' #' for more on explained variance). It is quite common to plot component 1 versus 2, then 1 versus 3
#' #' and finally 2 versus 3.
#' #' 
#' #' The argument \samp{show.labels} can be used to turn off the display of labels, only markers (dots)
#' #' will appear.
#' #' 
#' #' In \code{\link{plotScores}} you can specify alternative labels in \samp{labels}. By default, the
#' #' GID-tag is used for each genome. You can supply a vector of alternative labels. The labels may be
#' #' in any order, but the vector must be named by the GID-tags, i.e. each element in \samp{labels} must
#' #' have a name which is a valid GID-tag for some genome. This is necessary to ensure the alternative
#' #' labels are placed correctly in the score-space.
#' #' 
#' #' There is no alternative labelling of loading-plots, since the gene clusters lack a GID-tag-like system.
#' #' You can, however, change the gene cluster names by editing the column names of the pan-matrix directly
#' #' before you do the \code{\link{panpca}}.
#' #' 
#' #' You may color each label/marker individually. In \code{\link{plotScores}} you can again supply a vector
#' #' of colors, and name every element with a GID-tag to make certain they are used correctly. In
#' #' \code{\link{plotLoadings}} you can supply a vector of colors, but you must arrange them in proper
#' #' order yourself.
#' #' 
#' #' Additional arguments are passed on to \code{\link{text}} if \samp{show.labels=TRUE} and to
#' #' \code{\link{points}} if \samp{show.labels=FALSE}.
#' #' 
#' #' @author Lars Snipen and Kristian Hovde Liland.
#' #' 
#' #' @seealso \code{\link{panpca}}, \code{\link{plot.Panpca}}.
#' #' 
#' #' @examples 
#' #' # Loading a Panmat object in the micropan package
#' #' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' #' ppca.blast <- panpca(Mpneumoniae.blast.panmat)
#' #' 
#' #' # Plotting scores and loadings
#' #' plotScores(ppca.blast) # A score-plot
#' #' plotLoadings(ppca.blast) # A loading plot
#' #' 
#' #' # Plotting score with alternative labels and colors
#' #' data(list="Mpneumoniae.table",package="micropan")
#' #' labels <- Mpneumoniae.table$Strain
#' #' names(labels) <- Mpneumoniae.table$GID.tag
#' #' cols <- Mpneumoniae.table$Color
#' #' names(cols) <- Mpneumoniae.table$GID.tag
#' #' plotScores(ppca.blast,labels=labels,col=cols)
#' #' 
#' #' @export
#' plotScores <- function( pan.pca, x=1, y=2, show.labels=TRUE, labels=NULL, col="black", pch=16, ... ){
#'   Z <- pan.pca$Scores
#'   xr <- range( Z[,x] )
#'   yr <- range( Z[,y] )
#'   args <- list(...)
#'   ii <- match( c("xlab","ylab"), names( args ) )
#'   if( is.na(ii[1]) ){
#'     xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
#'   } else {
#'     xlab <- args[[ii[1]]]
#'   }
#'   if( is.na(ii[2]) ){
#'     ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
#'   } else {
#'     ylab <- args[[ii[2]]]
#'   }
#'   
#'   plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
#'   points( c(0,0), yr, type="l", col="gray" )
#'   
#'   pm.gid <- rownames( Z )
#'   if( length( col )>1 ){
#'     if( is.null( names( col ) ) ) stop( "Each element in col must be named by its GID.tag" )
#'     gid <- names( col )
#'     idx <- match( pm.gid, gid )
#'     if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
#'     cols <- col[idx]
#'     if( length( cols ) != length( pm.gid ) ) stop( "The number of elements in col does not correspond to the number of genomes in the pan-matrix" )
#'   } else {
#'     cols <- rep( col, length.out=length( pm.gid ) )
#'   }
#'   if( is.null( labels ) ){
#'     labs <- rownames( Z )
#'   } else {
#'     if( is.null( names( labels ) ) ) stop( "Each element in labels must be named by its GID.tag" )
#'     gid <- names( labels )
#'     idx <- match( pm.gid, gid )
#'     if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
#'     labs <- as.character( labels[idx] )
#'     if( length( labs ) != length( pm.gid ) ) stop( "The number of elements in labels does not correspond to the number of genomes in the pan-matrix" )
#'   }
#'   if( show.labels ){
#'     text( Z[,x], Z[,y], labs, col=cols, ... )
#'   } else {
#'     points( Z[,x], Z[,y], pch=pch, col=cols, ... )
#'   }
#' }
#' #' @rdname scores.loadings
#' #' @export
#' plotLoadings <- function( pan.pca, x=1, y=2, show.labels=TRUE, col="black", pch=16, ... ){
#'   L <- pan.pca$Loadings
#'   xr <- range( L[,x] )
#'   yr <- range( L[,y] )
#'   args <- list(...)
#'   ii <- match( c("xlab","ylab"), names( args ) )
#'   if( is.na(ii[1]) ){
#'     xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
#'   } else {
#'     xlab <- args[[ii[1]]]
#'   }
#'   if( is.na(ii[2]) ){
#'     ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
#'   } else {
#'     ylab <- args[[ii[2]]]
#'   }
#'   
#'   plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
#'   points( c(0,0), yr, type="l", col="gray" )
#'   
#'   if( show.labels ){
#'     text( L[,x], L[,y], rownames( L ), col=col, ... )
#'   } else {
#'     points( L[,x], L[,y], pch=pch, col=col, ... )
#'   }
#' }
