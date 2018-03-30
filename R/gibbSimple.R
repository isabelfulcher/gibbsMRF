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
#' @param burnin burnin
#' @param weights number of neighbors
#' @param start vector of starting values
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' N <- 100
#' x <- matrix(rbinom(N*N,1,0.1),N,N)
#' alpha <- c(.2,.5)
#' weights <- apply(x,1,sum)
#' output <- gibbSimple(alpha, x, N, 100, 10, weights, rbinom(N,1,runif(1,0,1)))
#'
#' @export
setGeneric(name = "gibbSimple", def = function(alpha, x, N, R, burnin, weights, start)
  standardGeneric("gibbSimple"))

#' @rdname gibbSimple
setMethod("gibbSimple", signature(alpha="vector", x="matrix", N="numeric", R="numeric", burnin="numeric", weights="vector", start="vector"),
          definition = function(alpha, x, N, R=500, burnin=0, weights, start) {

            stopifnot(dim(alpha) == 2)
            stopifnot(ncol(x) == nrow(x))

            iter <- R + burnin

            if (R==1){
              mrf <- run_gibbSimple1(x, weights, alpha[1], alpha[2], iter, N, start)
              mrfout <- mrf
            } else{
              mrf <- run_gibbSimple(x, weights, alpha[1], alpha[2], iter, N, start)
              mrfout <- mrf[,(burnin+1):iter]
            }

            return(mrfout)

          })
