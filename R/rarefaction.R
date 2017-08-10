#' @name rarefaction
#' @title Rarefaction curves for a pan-genome
#' 
#' @description Computes rarefaction curves for a number of random permutations of genomes.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param n.perm The number of random genome orderings to use. If \samp{n.perm=1} the fixed order of
#' the genomes in \samp{pan.matrix} is used.
#' 
#' @details A rarefaction curve is simply the cumulative number of unique gene clusters we observe as
#' more and more genomes are being considered. The shape of this curve will depend on the order of the
#' genomes. This function will typically compute rarefaction curves for a number of (\samp{n.perm})
#' orderings. By using a large number of permutations, and then averaging over the results, the effect
#' of any particular ordering is smoothed away.
#' 
#' The averaged curve illustrates how many new gene clusters we observe for each new genome. If this
#' levels out and becomes flat, it means we expect few, if any, new gene clusters by sequencing more
#' genomes. The function \code{\link{heaps}} can be used to estimate population openness based on this
#' principle.
#' 
#' @return This function returns a \code{Rarefac} object, which is a small extension to a matrix. The
#' generic functions \code{\link{plot.Rarefac}} and \code{\link{summary.Rarefac}}
#' are available for such objects.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{heaps}}, \code{\link{panMatrix}}, \code{\link{plot.Rarefac}},
#' \code{\link{summary.Rarefac}}.
#' 
#' @examples 
#' # Loading two Panmat objects in the micropan package 
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Rarefaction based on a BLAST clustering Panmat object
#' rarefac.blast <- rarefaction(Mpneumoniae.blast.panmat,n.perm=100)
#' plot(rarefac.blast)
#' 
#' # Rarefaction based on domain sequence clustering Panmat object
#' rarefac.domains <- rarefaction(Mpneumoniae.domain.panmat,n.perm=1000)
#' summary(rarefac.domains)
#' 
#' @export
rarefaction <- function( pan.matrix, n.perm=1 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  nmat <- matrix( 0, nrow=ng, ncol=n.perm )
  cm <- apply( pan.matrix, 2, cumsum )
  nmat[,1] <- rowSums( cm>0 )
  if( n.perm > 1 ){
    cat( "permuting:\n" )
    for( i in 2:n.perm ){
      cm <- apply( pan.matrix[sample( ng ),], 2, cumsum )
      nmat[,i] <- rowSums( cm>0 )
      cat( "." )
      if( (i/100)==round(i/100) ) cat( "\n" )
    }
    cat( "\n" )
  }
  rownames( nmat ) <- paste( 1:ng, "genomes" )
  colnames( nmat ) <- paste( "Permutation", 1:n.perm )
  class( nmat ) <- c( "Rarefac", "matrix" )
  return( nmat )
}



#' @rdname generic.Rarefac
#' @name plot.Rarefac
#' @title Plot and summary of \code{Rarefac} objects
#' 
#' @description Generic functions for \code{Rarefac} object.
#' 
#' @param x A \code{Rarefac} object, see below.
#' @param object A \code{Rarefac} object, see below.
#' @param type Type of plot, default is \samp{"b"}, giving markers with lines between.
#' @param pch Marker type, default is \samp{16}, a filled circle.
#' @param xlab Text for horizontal axis.
#' @param ylab Text for vertical axis.
#' @param \dots Optional graphical arguments.
#' 
#' @details A \code{Rarefac} object is a small (S3) extension to a matrix. The first column contains
#' the cumulative number of unique gene clusters found when considering 1,2,...,G genomes in a pan-matrix.
#' Thus, the \code{Rarefac} object is a matrix with G rows. Any additional columns will hold similar
#' numbers, but for random shufflings of the genome's ordering. A \code{Rarefac} object is typically
#' created by the function \code{\link{rarefaction}}.
#' 
#' The \code{\link{plot.Rarefac}} function will display the content of the \code{Rarefac} object as a plot
#' of the mean value in rows 1,2,...,G, where G is the total number of genomes in the study.
#' 
#' The \code{\link{summary.Rarefac}} function will display a text giving the same information as
#' \code{\link{plot.Rarefac}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{rarefaction}}, \code{\link{heaps}}.
#' 
#' @examples # See examples in the Help-file for rarefaction.
#' 
#' @export
plot.Rarefac <- function( x, type="b", pch=16, xlab="Genomes", ylab="Number of unique gene clusters", ... ){
  Rarefac <- x
  plot( 1:dim( Rarefac )[1], rowMeans( Rarefac ), type=type, pch=pch, xlab=xlab, ylab=ylab, ... )
}
#' @rdname generic.Rarefac
#' @export
summary.Rarefac <- function( object, ... ){
  cat( "For", 1, "genome we observe on average", round( mean( object[1,] ) ), "unique gene clusters\n" )
  for( i in 2:nrow( object ) ){
    cat( "For", i, "genomes we observe on average", round( mean( object[i,] ) ), "unique gene clusters\n" )
  }
}
