#' @name chao
#' @title The Chao lower bound estimate of pan-genome size
#' 
#' @description Computes the Chao lower bound estimated number of gene clusters in a pan-genome.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' 
#' @details The size of a pan-genome is the number of gene clusters in it, both those observed and those
#' not yet observed.
#' 
#' The input \samp{pan.matrix} is a a matrix with one row for each
#' genome and one column for each observed gene cluster in the pan-genome. See \code{\link{panMatrix}}
#' for how to construct this.
#' 
#' The number of observed gene clusters is simply the number of columns in \samp{pan.matrix}. The
#' number of gene clusters not yet observed is estimated by the Chao lower bound estimator (Chao, 1987).
#' This is based solely on the number of clusters observed in 1 and 2 genomes. It is a very simple and
#' conservative estimator, i.e. it is more likely to be too small than too large. 
#' 
#' @return The function returns an integer, the estimated pan-genome size. This includes both the number
#' of gene clusters observed so far, as well as the estimated number not yet seen.
#' 
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal
#' catchability. Biometrics, 43:783-791.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{binomixEstimate}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package
#' data(xmpl.panmat)
#' 
#' # Estimating the pan-genome size using the Chao estimator
#' chao.pansize <- chao(xmpl.panmat)
#' 
#' @export chao
#' 
chao <- function(pan.matrix){  
  y <- table(factor(colSums(pan.matrix > 0), levels = 1:nrow(pan.matrix)))
  if(y[2] == 0){
    stop( "Cannot compute Chao estimate since there are 0 gene clusters observed in 2 genomes!\n" )
  } else {
    pan.size <- round(sum(y) + y[1]^2/(2*y[2]))
    names(pan.size) <- NULL
    return(pan.size)
  }
}


#' @name heaps
#' @title Heaps law estimate
#' 
#' @description Estimating if a pan-genome is open or closed based on a Heaps law model.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param n.perm The number of random permutations of genome ordering.
#' 
#' @details An open pan-genome means there will always be new gene clusters observed as long as new genomes
#' are being sequenced. This may sound controversial, but in a pragmatic view, an open pan-genome indicates
#' that the number of new gene clusters to be observed in future genomes is \sQuote{large} (but not literally
#' infinite). Opposite, a closed pan-genome indicates we are approaching the end of new gene clusters. 
#' 
#' This function is based on a Heaps law approach suggested by Tettelin et al (2008). The Heaps law model
#' is fitted to the number of new gene clusters observed when genomes are ordered in a random way. The model
#' has two parameters, an intercept and a decay parameter called \samp{alpha}. If \samp{alpha>1.0} the
#' pan-genome is closed, if \samp{alpha<1.0} it is open.
#' 
#' The number of permutations, \samp{n.perm}, should be as large as possible, limited by computation time.
#' The default value of 100 is certainly a minimum.
#' 
#' Word of caution: The Heaps law assumes independent sampling. If some of the genomes in the data set
#' form distinct sub-groups in the population, this may affect the results of this analysis severely.
#' 
#' @return A vector of two estimated parameters: The \samp{Intercept} and the decay parameter \samp{alpha}.
#' If \samp{alpha<1.0} the pan-genome is open, if \samp{alpha>1.0} it is closed.
#' 
#' @references Tettelin, H., Riley, D., Cattuto, C., Medini, D. (2008). Comparative genomics: the
#' bacterial pan-genome. Current Opinions in Microbiology, 12:472-477.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{binomixEstimate}}, \code{\link{chao}}, \code{\link{rarefaction}}.
#' 
#' @examples 
#' # Loading a pan-matrix in this package 
#' data(xmpl.panmat)
#' 
#' # Estimating population openness
#' h.est <- heaps(xmpl.panmat, n.perm = 500)
#' # If alpha < 1 it indicates an open pan-genome
#' 
#' @importFrom stats optim
#' 
#' @export
heaps <- function(pan.matrix, n.perm = 100){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  ng <- dim(pan.matrix)[1]
  nmat <- matrix(0, nrow = nrow(pan.matrix) - 1, ncol = n.perm)
  for(i in 1:n.perm){
    cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)
    nmat[,i] <- rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])
    cat(i, "/", n.perm, "\r")
  }
  x <- rep((2:nrow(pan.matrix)), times = n.perm)
  y <- as.numeric(nmat)
  p0 <- c(mean(y[which(x == 2)] ), 1)
  fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))
  p.hat <- fit$par
  names(p.hat) <- c("Intercept", "alpha")
  return(p.hat)
}

objectFun <- function(p, x, y){
  y.hat <- p[1] * x^(-p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}
