#' @include gibbsMRF.R
NULL

#' Execute the gibbs sampler and save sufficient statistics for a simple case
#'
#' \code{gibbSimplestat} saves sufficient statistics
#' for a Markov Random Field model. This is done
#' iteratively in \code{Rcpp} to minimize the computing time.
#'
#' @param x adjacency matrix for network
#' @param alpha vector of alpha parameters
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
#' R <- 100
#' output <- gibbSimplestat(x, alpha, N, R, weights)
#'
#' @export
setGeneric(name = "gibbSimplestat", def = function(x, alpha, N, R, weights)
  standardGeneric("gibbSimplestat"))

#' @rdname gibbSimplestat
setMethod("gibbSimplestat", signature(x="matrix", alpha="vector", N="numeric", R="numeric", weights="vector"),
          definition = function(x, alpha, N, R, weights) {

            stopifnot(dim(alpha) == 2)
            stopifnot(ncol(x) == nrow(x))

            start <- rbinom(N,1,runif(1,0,1)) #initialize starting vector
            statout <- run_gibbSimplestat(x, weights, alpha[1], alpha[2], R, N, start)

            return(statout)

          })
