#' @include gibbsMRF.R
NULL

#' Execute the gibbs sampler for a simple case
#'
#' \code{gibbSimple} creates multiple Markov Random Field realizations
#' for one covariate model. This is done iteratively in \code{Rcpp}
#' to minimize the computing time.
#'
#' @param x adjacency matrix for network
#' @param alpha parameters
#' @param N number of nodes
#' @param R number of iterations
#' @param weights number of neighbors
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' N <- 100
#' x <- matrix(rbinom(N*N,1,0.1),N,N)
#' alpha <- c(.2,.5)
#' weights <- apply(x,1,sum)
#' R <- 1000
#' output <- gibbSimple(x, alpha, N, R, weights)
#'
#' @export
setGeneric(name = "gibbSimple", def = function(x, alpha, N, R, weights)
  standardGeneric("gibbSimple"))

#' @rdname gibbSimple
setMethod("gibbSimple", signature(x="matrix", alpha="vector", N="numeric", R="numeric", weights="vector"),
          definition = function(x, alpha, N, R, weights) {

            stopifnot(dim(alpha) == 2)
            stopifnot(ncol(x) == nrow(x))

            start <- rbinom(N,1,runif(1,0,1)) #initialize starting vector
            mrfout <- run_gibbSimple(x, weights, alpha[1], alpha[2], R, N, start)

            return(mrfout)

          })
