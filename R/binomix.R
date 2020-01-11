#' @name binomixEstimate
#' @aliases binomixEstimate
#' 
#' @title Binomial mixture model estimates
#' 
#' @description Fits binomial mixture models to the data given as a pan-matrix. From the fitted models
#' both estimates of pan-genome size and core-genome size are available.
#' 
#' @param pan.matrix A pan-matrix, see \code{\link{panMatrix}} for details.
#' @param K.range The range of model complexities to explore. The vector of integers specify the number
#' of binomial densities to combine in the mixture models.
#' @param core.detect.prob The detection probability of core genes. This should almost always be 1.0,
#' since a core gene is by definition always present in all genomes, but can be set fractionally smaller.
#' @param verbose Logical indicating if textual output should be given to monitor the progress of the
#' computations.
#' 
#' @details  A binomial mixture model can be used to describe the distribution of gene clusters across
#' genomes in a pan-genome. The idea and the details of the computations are given in Hogg et al (2007),
#' Snipen et al (2009) and Snipen & Ussery (2012).
#' 
#' Central to the concept is the idea that every gene has a detection probability, i.e. a probability of
#' being present in a genome. Genes who are always present in all genomes are called core genes, and these
#' should have a detection probability of 1.0. Other genes are only present in a subset of the genomes, and
#' these have smaller detection probabilities. Some genes are only present in one single genome, denoted
#' ORFan genes, and an unknown number of genes have yet to be observed. If the number of genomes investigated
#' is large these latter must have a very small detection probability. 
#' 
#' A binomial mixture model with \samp{K} components estimates \samp{K} detection probabilities from the
#' data. The more components you choose, the better you can fit the (present) data, at the cost of less
#' precision in the estimates due to less degrees of freedom. \code{\link{binomixEstimate}} allows you to
#' fit several models, and the input \samp{K.range} specifies which values of \samp{K} to try out. There no
#' real point using \samp{K} less than 3, and the default is \samp{K.range=3:5}. In general, the more genomes
#' you have the larger you can choose \samp{K} without overfitting.  Computations will be slower for larger
#' values of \samp{K}. In order to choose the optimal value for \samp{K}, \code{\link{binomixEstimate}}
#' computes the BIC-criterion, see below.
#' 
#' As the number of genomes grow, we tend to observe an increasing number of gene clusters. Once a
#' \samp{K}-component binomial mixture has been fitted, we can estimate the number of gene clusters not yet
#' observed, and thereby the pan-genome size. Also, as the number of genomes grows we tend to observe fewer
#' core genes. The fitted binomial mixture model also gives an estimate of the final number of core gene
#' clusters, i.e. those still left after having observed \sQuote{infinite} many genomes.
#' 
#' The detection probability of core genes should be 1.0, but can at times be set fractionally smaller.
#' This means you accept that even core genes are not always detected in every genome, e.g. they may be
#' there, but your gene prediction has missed them. Notice that setting the \samp{core.detect.prob} to less
#' than 1.0 may affect the core gene size estimate dramatically.
#' 
#' @return \code{\link{binomixEstimate}} returns a \code{list} with two components, the \samp{BIC.tbl}
#' and \samp{Mix.tbl}.
#' 
#' The \samp{BIC.tbl} is a \code{tibble} listing, in each row, the results for each number of components
#' used, given by the input \samp{K.range}. The column \samp{Core.size} is the estimated number of
#' core gene families, the column \samp{Pan.size} is the estimated pan-genome size. The column
#' \samp{BIC} is the Bayesian Information Criterion (Schwarz, 1978) that should be used to choose the
#' optimal component number (\samp{K}). The number of components where \samp{BIC} is minimized is the
#' optimum. If minimum \samp{BIC} is reached for the largest \samp{K} value you should extend the
#' \samp{K.range} to larger values and re-fit. The function will issue
#' a \code{warning} to remind you of this.
#' 
#' The \samp{Mix.tbl} is a \code{tibble} with estimates from the mixture models. The column \samp{Component}
#' indicates the model, i.e. all rows where \samp{Component} has the same value are from the same model.
#' There will be 3 rows for 3-component model, 4 rows for 4-component, etc. The column \samp{Detection.prob}
#' contain the estimated detection probabilities for each component of the mixture models. A 
#' \samp{Mixing.proportion} is the proportion of the gene clusters having the corresponding \samp{Detection.prob},
#' i.e. if core genes have \samp{Detection.prob} 1.0, the corresponding \samp{Mixing.proportion} (same row)
#' indicates how large fraction of the gene families are core genes.
#' 
#' @references
#' Hogg, J.S., Hu, F.Z, Janto, B., Boissy, R., Hayes, J., Keefe, R., Post, J.C., Ehrlich, G.D. (2007).
#' Characterization and modeling of the Haemophilus influenzae core- and supra-genomes based on the
#' complete genomic sequences of Rd and 12 clinical nontypeable strains. Genome Biology, 8:R103.
#' 
#' Snipen, L., Almoy, T., Ussery, D.W. (2009). Microbial comparative pan-genomics using binomial
#' mixture models. BMC Genomics, 10:385.
#' 
#' Snipen, L., Ussery, D.W. (2012). A domain sequence approach to pangenomics: Applications to
#' Escherichia coli. F1000 Research, 1:19.
#' 
#' Schwarz, G. (1978). Estimating the Dimension of a Model. The Annals of Statistics, 6(2):461-464.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{chao}}.
#' 
#' @examples
#' # Loading an example pan-matrix
#' data(xmpl.panmat)
#' 
#' # Estimating binomial mixture models
#' binmix.lst <- binomixEstimate(xmpl.panmat, K.range = 3:8)
#' print(binmix.lst$BIC.tbl) # minimum BIC at 3 components
#' 
#' # The pan-genome gene distribution as a pie-chart
#' ncomp <- 3
#' binmix.lst$Mix.tbl %>% 
#'   filter(Components == ncomp) %>% 
#'   ggplot() +
#'   geom_col(aes(x = "", y = Mixing.proportion, fill = Detection.prob)) +
#'   coord_polar(theta = "y") +
#'   labs(x = "", y = "", title = "Pan-genome gene distribution") +
#'   scale_fill_gradientn(colors = c("pink", "orange", "green", "cyan", "blue"))
#'   
#' # The distribution in an average genome
#' binmix.lst$Mix.tbl %>% 
#'   filter(Components == ncomp) %>% 
#'   mutate(Single = Mixing.proportion * Detection.prob) %>%
#'   ggplot() +
#'   geom_col(aes(x = "", y = Single, fill = Detection.prob)) +
#'   coord_polar(theta = "y") +
#'   labs(x = "", y = "", title = "Average genome gene distribution") +
#'   scale_fill_gradientn(colors = c("pink", "orange", "green", "cyan", "blue"))
#' 
#' @importFrom stringr str_c
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot
#' 
#' @export binomixEstimate
#' 
binomixEstimate <- function(pan.matrix, K.range = 3:5, core.detect.prob = 1.0, verbose = TRUE){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  y <- table(factor(colSums(pan.matrix), levels = 1:nrow(pan.matrix)))
  bic.tbl <- tibble(Components = K.range,
                    Core.size = rep(0, length(K.range)),
                    Pan.size = rep(0, length(K.range)),
                    BIC = rep(0, length(K.range)))
  mix.tbl <- NULL
  for(i in 1:length(K.range)){
    if(verbose) cat("binomixEstimate: Fitting", K.range[i], "component model...\n")
    lst <- binomixMachine(y, K.range[i], core.detect.prob)
    bic.tbl[i,-1] <- lst[[1]]
    mix.tbl <- bind_rows(mix.tbl, lst[[2]])
  }
  if(bic.tbl[length(K.range),3] == min(bic.tbl[,3])) warning("Minimum BIC at maximum K, increase upper limit of K.range")
  return(list(BIC.tbl = bic.tbl, Mix.tbl = mix.tbl))
}


#' @importFrom stats constrOptim as.dendrogram dendrapply dist is.leaf optim prcomp sd
#' @importFrom tibble tibble
binomixMachine <- function(y, K, core.detect.prob = 1.0){
  n <- sum(y)
  G <- length(y)
  ctr <- list(maxit = 300, reltol = 1e-6)
  np <- K - 1
    
  pmix0 <- rep(1, np)/K             # flat mixture proportions
  pdet0 <- (1:np)/(np+1)            # "all" possible detection probabilities
  p.initial <- c(pmix0, pdet0)      # initial values for parameters
  # the inequality constraints...
  A <- rbind(c( rep(1, np), rep(0, np)), c(rep(-1, np), rep(0, np)), diag(np+np), -1*diag(np+np))
  b <- c(0, -1, rep(0, np+np), rep(-1, np+np))
  
  # The estimation, maximizing the negative truncated log-likelihood function
  est <- constrOptim(theta = p.initial, f = negTruncLogLike, grad = NULL, method = "Nelder-Mead", control = ctr, ui = A, ci = b,
                     y = y, core.p = core.detect.prob)
  
  estimates <- numeric(3)
  names(estimates) <- c("Core.size", "Pan.size", "BIC")
  estimates[3] <- 2*est$value + log(n)*(np+K)                       # the BIC-criterion
  p.mix <- c(1 - sum(est$par[1:np]), est$par[1:np])                 # the mixing proportions
  p.det <- c(core.detect.prob, est$par[(np+1):length( est$par )])   # the detection probabilities
  ixx <- order(p.det)
  p.det <- p.det[ixx]
  p.mix <- p.mix[ixx]
    
  theta_0 <- choose(G, 0) * sum(p.mix * (1-p.det)^G)
  y_0 <- n * theta_0/(1-theta_0)
  estimates[2] <- n + round(y_0)
  ixx <- which(p.det >= core.detect.prob)
  estimates[1] <- round(estimates[2] * sum(p.mix[ixx]))
  mix.tbl <- tibble(Components = rep(K, length(p.det)),
                    Detection.prob = p.det,
                    Mixing.proportion = p.mix)
  return(list(estimates, mix.tbl))
}

negTruncLogLike <- function(p, y, core.p){
  np <- length(p)/2
  p.det <- c(core.p, p[(np+1):length(p)])
  p.mix <- c(1-sum(p[1:np]), p[1:np])
  G <- length(y)
  K <- length(p.mix)
  n <- sum(y)
    
  theta_0 <- choose(G, 0) * sum(p.mix * (1-p.det)^G)
  L <- -n * log(1 - theta_0)
  for(g in 1:G){
    theta_g <- choose(G, g) * sum(p.mix * p.det^g * (1-p.det)^(G-g))
    L <- L + y[g] * log(theta_g)
  }
  return(-L)
}









#' #' @rdname generic.Binomix
#' #' @name plot.Binomix
#' #' 
#' #' @title Plot and summary of \code{Binomix} objects
#' #' 
#' #' @description Generic functions for \code{Binomix} objects.
#' #' 
#' #' @param x A \code{Binomix} object, see below.
#' #' @param object A \code{Binomix} object, see below.
#' #' @param type Type of plot, default is \samp{type="pan"} which means the pie chart shows distribution
#' #' over the entire pan-genome. The alternative is \samp{type="single"} which means the pie chart will
#' #' show the distribution within a single (average) genome.
#' #' @param cex Plot symbol scaling.
#' #' @param ncomp Which model to display. You can override the display of the optimal (minimum BIC) model
#' #' by specifying the number of components here, e.g. \samp{ncomp=5} will always display the model with
#' #' \samp{5} components regardless of its BIC value.
#' #' @param show.bar Logical indicating if a colorbar should be displayed next to the pie.
#' #' @param \dots Optional graphical arguments.
#' #' 
#' #' @details A \code{Binomix} object contains a series of fitted binomial mixture models. It is a small
#' #' (S3) extension to a \code{list}, having two components. These are named \samp{BIC.table} and
#' #' \samp{Mix.list}, see \code{\link{binomixEstimate}} for more details.
#' #' 
#' #' The \code{\link{plot.Binomix}} function will display a \code{Binomix} object as a pie chart. Only
#' #' the model with the smallest BIC-criterion value is displayed. The BIC-criterion is used to rank the
#' #' various fitted models, and minimum BIC is an objective criterion for finding the best model complexity.
#' #' Each sector of the pie chart is a component, the color of the sector indicates its detection probability
#' #' and the size of the sector its mixing proportion. This pie chart illustrates how gene clusters are
#' #' distributed within the pan-genome. Sectors of (dark) blue color are highly conserved gene clusters
#' #' (core genes), sectors of greenish colors are medium conserved clusters (shell genes) and sectors of
#' #' orange/pink colors are non-conserved clusters (cloud genes).
#' #' 
#' #' The \code{\link{summary.Binomix}} function will print the estimated core size and pan-genome size for
#' #' the optimal component model.
#' #' 
#' #' @author Lars Snipen and Kristian Hovde Liland.
#' #' 
#' #' @seealso \code{\link{binomixEstimate}}.
#' #' 
#' #' @examples # See examples in the Help-file for binomixEstimate.
#' #' 
#' #' @importFrom graphics axis barplot box layout par pie plot points rect text
#' #' 
#' #' @export
#' plot.Binomix <- function(x, type = "pan", cex = 2, ncomp = NA, show.bar = TRUE, ...){
#'   Binomix <- x
#'   if(is.na(ncomp)){
#'     ncomp <- which(Binomix$BIC.tbl[,3] == min(Binomix$BIC.tbl[,3]))[1]
#'   } else {
#'     ncomp <- which(as.numeric(str_remove(rownames(Binomix$BIC.tbl), " components")) == ncomp)
#'     if(length(ncomp) != 1) stop("Specified value of ncomp has not been fitted for this model")
#'   }
#'   cpar <- par()$mar
#'   par(mar=c(2,2,2,2))
#'   typ <- grep(type, c("pan", "single"))
#'   fit <- Binomix$Mix.lst[[ncomp]]
#'   dprob <- as.character(round(fit[1,]*1000)/1000)
#'   if(show.bar){
#'     layout(matrix(c(1,1,1,1,1,1,2), nrow=1))
#'     if(length(typ) == 0){
#'       stop("Unknown type specified")
#'     } else if(typ == 1){
#'       pie(fit[2,], clockwise = T, col = panColor(fit[1,]), labels = dprob, radius = 1.0, cex = cex, ...)
#'     } else {
#'       eg <- fit[2,]*fit[1,]
#'       pie(eg/sum(eg), clockwise = T, col = panColor(fit[1,]), labels = dprob, radius = 1.0, cex = cex, ...)
#'     }
#'     par(mar = c(1,1,1,3))
#'     p <- (0:100)
#'     plot(rep(0, 101), p, cex = 0, col = "white", xlim = c(0,1), ylim = c(0,100), 
#'           xaxt = "n", yaxt = "n", xlab= "", ylab = "")
#'     box(lwd = 3, col = "white")
#'     cols <- panColor(p/100)
#'     xl <- rep(0, 100)
#'     xr <- rep(1, 100)
#'     yb <- (0:100)
#'     yt <- 1:101
#'     rect(xl, yb, xr, yt, col = cols, border = cols)
#'     axis(side = 4, at = seq(0, 100, 10), labels = as.character(seq(0, 100, 10)/100))
#'   } else {
#'     if(length(typ) == 0){
#'       stop("Unknown type specified")
#'     } else if(typ == 1){
#'       pie(fit[2,], clockwise = T, col = panColor(fit[1,]), labels = dprob, radius = 1.0, cex = cex, ...)
#'     } else {
#'       eg <- fit[2,]*fit[1,]
#'       pie(eg/sum(eg), clockwise = T, col = panColor(fit[1,]), labels = dprob, radius = 1.0, cex = cex, ...)
#'     }
#'   }
#'   par(mar = cpar)
#' }
#' #' @rdname generic.Binomix
#' #' @export
#' summary.Binomix <- function(object, ...){
#'   ncomp <- which(object$BIC.tbl[,3] == min(object$BIC.tbl[,3]))[1]
#'   cat("Minimum BIC model at", rownames(object$BIC.tbl)[ncomp], "\nFor this model:\n")
#'   cat("Estimated core size:", object$BIC.tbl[ncomp,1], "clusters\n")
#'   cat("Estimated pangenome size:", object$BIC.tbl[ncomp,2], "clusters\n")
#' }
#' 
#' plotBar <- function(){
#'   cpar <- par()$mar
#'   par(mar = c(1,1,1,3))
#'   p <- (0:100)
#'   plot(rep(0, 101), p, cex = 0, col = "white", xlim = c(0,1), ylim = c(0,100), 
#'         xaxt = "n", yaxt = "n", xlab = "", ylab = "" )
#'   box(lwd = 3, col = "white")
#'   cols <- panColor(p/100)
#'   xl <- rep(0, 100)
#'   xr <- rep(1, 100)
#'   yb <- (0:100)
#'   yt <- 1:101
#'   rect(xl, yb, xr, yt, col = cols, border = cols)
#'   axis(side = 4, at = seq(0, 100, 10), labels = as.character(seq(0, 100, 10)/100))
#' }
#' 
#' #' @importFrom grDevices colorRampPalette
#' panColor <- function(p.vec){
#'   level <- pretty(c(0, 1), 100)
#'   nlevel <- length(level)
#'   crp <- colorRampPalette(c("pink", "orange", "green", "cyan", "blue"))(nlevel)
#'   return(crp[1+round(100*p.vec)])
#' }
#' 


