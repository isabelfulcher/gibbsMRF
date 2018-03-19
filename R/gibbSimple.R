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
#' @param R number of iterations
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' N <- 100
#' x <- matrix(rbinom(N*N,1,0.1),N,N)
#' alpha <- c(.2,.5)
#' R <- 1000
#' output <- gibbSimple(x, alpha, R)
#'
#' @export
setGeneric(name = "gibbSimple", def = function(x, alpha, R)
  standardGeneric("gibbSimple"))

#' @rdname gibbSimple
setMethod("gibbSimple", signature(x="matrix", alpha="vector", R="numeric"),
          definition = function(x, alpha, R=1000) {

            stopifnot(dim(alpha) == 2)
            stopifnot(ncol(x) == nrow(x))

            N <- nrow(x) #number of nodes
            weights <- apply(x,1,sum) #neighbor weights
            start <- rbinom(N,1,runif(1,0,1)) #initialize starting vector

            mrfout <- run_gibbSimple(x, weights, alpha[1], alpha[2], R, N, start)

            return(mrfout)

          })
