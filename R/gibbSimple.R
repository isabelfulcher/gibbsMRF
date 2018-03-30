#' @include gibbsMRF.R
NULL

#' Execute the gibbs sampler for a simple case
#'
#' \code{gibbSimple} creates multiple Markov Random Field realizations
#' for one covariate model. This is done iteratively in \code{Rcpp}
#' to minimize the computing time.
#'
#' @param alpha parameters
#' @param x adjacency matrix for network
#' @param N number of nodes
#' @param R number of output iterations
#' @param weights number of neighbors
#' @param burnin cutoff value for burnin
#' @param thin thinning value for chain
#' @param start vector of starting values
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' N <- 100
#' x <- matrix(rbinom(N*N,1,0.1),N,N)
#' alpha <- c(.2,.5)
#' weights <- apply(x,1,sum)
#' output <- gibbSimple(alpha, x, N, 100, weights, 10, 2, rbinom(N,1,runif(1,0,1)))
#'
#' @export
setGeneric(name = "gibbSimple", def = function(alpha, x, N, R, weights, burnin, thin, start)
  standardGeneric("gibbSimple"))

#' @rdname gibbSimple
setMethod("gibbSimple", signature(alpha="vector", x="matrix", N="numeric", R="numeric", weights="vector", burnin="numeric", thin="numeric", start="vector"),
          definition = function(alpha, x, N, R=500, weights, burnin=0, thin=1,start) {

            stopifnot(dim(alpha) == 2)
            stopifnot(ncol(x) == nrow(x))

            iter <- R*thin + burnin
            mrf <- run_gibbSimple(x, weights, alpha[1], alpha[2], iter, N, start)
            if (iter==1){ mrfout <- mrf} else { mrfout <- mrf[,(burnin+1):iter][,seq(1,R*thin,by=thin)]}

            return(mrfout)

          })
