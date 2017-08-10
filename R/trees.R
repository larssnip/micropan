#' @name panTree
#' @title Constructing pan-genome trees
#' 
#' @description Creates a pan-genome tree based on a pan-matrix and a distance function.
#' 
#' @param  pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}}.
#' @param dist.FUN A valid distance function, see below.
#' @param nboot Number of bootstrap samples.
#' @param linkage The linkage function, see below.
#' @param \dots Additional parameters passed on to the specified distance function, see Details below.
#' 
#' @details A pan-genome tree is a graphical display of the genomes in a pan-genome study, based on
#' some pan-matrix (Snipen & Ussery, 2010). \code{\link{panTree}} is a constructor that computes a
#' \code{Pantree} object, use \code{\link{plot.Pantree}} to actually plot the tree.
#' 
#' The parameter \samp{dist.FUN} must be a function that takes as input a numerical matrix (\code{Panmat}
#' object) and returns a \code{\link{dist}} object. See \code{\link{distManhattan}} or
#' \code{\link{distJaccard}} for examples of such functions. Any additional arguments (\samp{...}) are
#' passed on to this function.
#' 
#' If you want to have bootstrap-values in the tree, set \samp{nboot} to some appropriate number (e.g.
#' \samp{nboot=100}).
#' 
#' The tree is created by \code{\link{hclust}} (hierarchical clustering) using the \samp{average}
#' linkage function, which is according to Snipen & Ussery, 2010. You may specify alternatives by the
#' parameter \samp{linkage}, see \code{\link{hclust}} for details.
#' 
#' @return  This function returns a \code{Pantree} object, which is a small (S3) extension to a
#' \code{\link{list}} with 4 components. These components are named \samp{Htree}, \samp{Nboot},
#' \samp{Nbranch} and \samp{Dist.FUN}.
#' 
#' \samp{Htree} is a \code{\link{hclust}} object. This is the actual tree.
#' \samp{Nboot} is the number of bootstrap samples.
#' \samp{Nbranch} is a vector listing the number of times each split/clade in the tree was observed
#' in the bootstrap procedure.
#' \samp{Dist.FUN} is the name of the distance function used to construct the tree.
#' 
#' @references Snipen, L., Ussery, D.W. (2010). Standard operating procedure for computing pangenome
#' trees. Standards in Genomic Sciences, 2:135-141.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distManhattan}}, \code{\link{distJaccard}},
#' \code{\link{plot.Pantree}}.
#' 
#' @examples 
#' # Loading a Panmat object, constructing a tree and plotting it 
#' data(list="Mpneumoniae.blast.panmat",package="micropan")
#' my.tree <- panTree(Mpneumoniae.blast.panmat)
#' plot(my.tree)
#' 
#' # Computing some weights to be used in the distManhattan
#' # function below...
#' w <- geneWeights(Mpneumoniae.blast.panmat,type="shell")
#' # Creating another tree with scaled and weighted distances and bootstrap values
#' my.tree <- panTree(Mpneumoniae.blast.panmat, scale=0.1, weights=w)
#' 
#' # ...and plotting with alternative labels and colors from Mpneumoniae.table
#' data(list="Mpneumoniae.table",package="micropan")
#' labels <- Mpneumoniae.table$Strain
#' names(labels) <- Mpneumoniae.table$GID.tag
#' cols <- Mpneumoniae.table$Color
#' names(cols) <- Mpneumoniae.table$GID.tag
#' plot(my.tree, leaf.lab=labels, col=cols,cex=0.8, xlab="Shell-weighted Manhattan distances")
#' 
#' @export
panTree <- function( pan.matrix, dist.FUN=distManhattan, nboot=0, linkage="average", ... ){
  P <- ncol( pan.matrix )
  distFUN <- match.fun( dist.FUN )
  htree <- hclust( distFUN( pan.matrix, ... ), method=linkage )
  
  nbranch=NULL
  if( nboot > 0 ){
    cat( "bootstrapping" )
    signatur <- clusterSignature( htree$merge )
    nbranch <- rep( 0, length( signatur ) )
    names( nbranch ) <- signatur
    for( i in 1:nboot ){
      cat( "." )
      ht <- hclust( distFUN( pan.matrix[,sample( (1:P), P, replace=T )], ... ), method=linkage )
      nbranch <- nbranch + as.numeric( signatur %in% clusterSignature( ht$merge ) )
    }
    cat( "\n" )
  }
  pantree <- list( Htree=htree, Nboot=nboot, Nbranch=nbranch, Dist.FUN=as.character( substitute( dist.FUN ) ) )
  class( pantree ) <- c( "Pantree", "list" )
  return( pantree )
}


#' @rdname generic.Pantree
#' @name plot.Pantree
#' @title Plot and summary of \code{Pantree} objects
#' 
#' @description Generic functions for \code{Pantree} objects.
#' 
#' @param x A \code{Pantree} object, see below.
#' @param object A \code{Pantree} object, see below.
#' @param leaf.lab Alternative labels for the leaves, see below.
#' @param col Color of the leaf labels, see below.
#' @param xlab Text for the x-axis (distance-axis) of the plotted tree.
#' @param main Title above the plotted tree.
#' @param cex Scaling of the leaf labels of the plotted tree.
#' @param show.boot Logical to turn off plotting of bootstrap values.
#' @param \dots Additional arguments, see below.
#' 
#' @details  A \code{Pantree} object is created by \code{\link{panTree}} and contains information to
#' display a pan-genome tree. The \code{\link{plot.Pantree}} function will display the tree as a
#' \code{\link{dendrogram}} object.
#' 
#' The argument \samp{leaf.lab} can be used to give alternative labels, the GID-tags are used by
#' default. \samp{leaf.lab} must be a vector of labels, one for each genome in the \code{Pantree}. The
#' labels may be in any order, but the vector must be named by the GID-tags, i.e. each element in
#' \samp{leaf.lab} must have a name which is a valid GID-tag for some genome. This is necessary to ensure
#' the alternative labels are placed correctly in the tree.
#' 
#' The argument \samp{col} specifies the color(s) of the leaf labels in the tree. It can either be a single
#' color or a vector of colors, one for each leaf label (genome). Again, the colors may be in any order,
#' but the vector must be named by the GID-tags, i.e. each element in \samp{col} must have a name which is
#' a valid GID-tag for some genome.
#' 
#' The argument \samp{cex} scales the leaf label font size.
#' 
#' The argument \samp{show.boot} can be used to turn off the display of bootstrap values. Note that if
#' the tree was constructed without bootstrapping, no bootstrap values are available, and this argument has
#' no effect.
#' 
#' Any additional arguments are passed on to the \code{\link{plot.dendrogram}} function.
#' 
#' \code{\link{summary.Pantree}} prints a short text describing
#' the \code{Pantree} object.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @note Using \samp{nodePar} to manipulate the \code{\link{dendrogram}} object will have no effect
#' on the leaf nodes here since these are set by the \code{\link{dendrapply}} function. The tree is
#' always displayed horizontal, to align the labels in the right margin for easy reading.
#' 
#' @seealso \code{\link{panTree}}.
#' 
#' @examples # See examples in the Help-file for panTree.
#' 
#' @export
plot.Pantree <- function( x, leaf.lab=NULL, col="black", xlab="", main="", cex=1, show.boot=TRUE, ... ){
  if( is.null( leaf.lab ) ){
    labs <- x$Htree$labels
  } else {
    if( is.null( names( leaf.lab ) ) ) stop( "Each element in leaf.lab must be named by its GID.tag" )
    pm.gid <- x$Htree$labels
    gid <- names( leaf.lab )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the tree" )
    labs <- as.character( leaf.lab[idx] )
    if( length( labs ) != length( pm.gid ) ) stop( "The number of elements in leaf.lab does not correspond to the number of genomes in the tree" )
  }
  # colors...
  if( length( col )==1 ){
    col.lab <- rep( col, length.out=length( labs ) )
  } else {
    if( is.null( names( col ) ) ) stop( "Each element in col must be named by its GID.tag" )
    pm.gid <- x$Htree$labels
    gid <- names( col )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the tree" )
    col.lab <- col[idx]
    if( length( col ) != length( pm.gid ) ) stop( "The number of elements in col does not correspond to the number of genomes in the tree" )
  }
  x$Htree$labels <- labs
  dendro <- dendrapply( as.dendrogram( x$Htree ), setLeafAttributes, labs, col.lab, cex )
  
  cpar <- par()$mar
  maxnc <- max( nchar( labs ) )
  mar4 <- 2 + cex*(maxnc*0.8*par("cin")[1]*2.54)
  par( mar=c(5,1,1,mar4) )
  plot( dendro, horiz=T, xlab=xlab, main=main, ... )
  if( (x$Nboot > 0) & show.boot ){
    lab.nod <- as.character( round( 100*x$Nbranch/x$Nboot )/100 )
    bp <- branchPos( x$Htree$merge, x$Htree$order )
    for( i in 1:length( x$Htree$height ) ){
      if( !(x$Htree$merge[i,1]<0 & x$Htree$merge[i,2]<0) ){
        text( x$Htree$height[i], bp[i], lab.nod[i], cex=0.75, col="red4", pos=4, offset=0.1 )
      }
    }
  }
  par( mar=cpar )  
}
#' @rdname generic.Pantree
#' @export
summary.Pantree <- function( object, ... ){
  labs <- object$Htree$labels
  cat( "Pangenome tree for ", length( labs ), " genomes (", paste( labs[1:3], collapse="," ), "...)\n", sep="" )
  cat( "Distances computed by", object$Dist.FUN, "and containing", object$Nboot, "bootstrap samples\n")
}

setLeafAttributes <- function( node, lab.names, lab.col, lab.cex ){
  if( is.leaf( node ) ){
    attr( node, "nodePar" ) <- list( lab.cex=lab.cex, pch=NA, lab.col=lab.col[which( lab.names == attr( node, "label" ) )] )
  }
  return( node )
}

clusterSignature <- function( mergeMatrix ){
  N <- nrow( mergeMatrix )
  signature <- character( N )
  for( i in 1:N ){
    if( mergeMatrix[i,1] < 0 ){
      left <- as.character( -1*mergeMatrix[i,1] )
    } else {
      left <- gsub( ";", ",", signature[mergeMatrix[i,1]] )
    }
    if( mergeMatrix[i,2] < 0 ){
      right <- as.character( -1*mergeMatrix[i,2] )
    } else {
      right <- gsub( ";", ",", signature[mergeMatrix[i,2]] )
    }
    left <- paste( sort( unlist( strsplit( left, split="," ) ) ), collapse="," )
    right <- paste( sort( unlist( strsplit( right, split="," ) ) ), collapse="," )
    zig <- sort( c( left, right ) )
    signature[i] <- paste( zig[1], zig[2], sep=";" )
  }
  return( signature )
}

branchPos <- function( mergeMatrix, ordering ){
  N <- nrow( mergeMatrix )
  branchp <- numeric( N )
  for( i in 1:N ){
    if( mergeMatrix[i,1] < 0 ){
      left <- which( ordering == -1*mergeMatrix[i,1] )
    } else {
      left <- branchp[mergeMatrix[i,1]]
    }
    if( mergeMatrix[i,2] < 0 ){
      right <- which( ordering == -1*mergeMatrix[i,2] )
    } else {
      right <- branchp[mergeMatrix[i,2]]
    }
    branchp[i] <- mean( c( left, right ) )
  }
  return( branchp )
}
