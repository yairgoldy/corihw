#' toy_example: Creates a toy example data of multiple correlated tests
#'
#' Create \emph{m} two-sample t-tests with correlations. Let \emph{mu} be the mean under
#' the alternative and 0 the mean under the null. The proportion of alternatives is \emph{prob}.
#' The first sample of \emph{n}  \emph{m}-dimendional observations is drawn from a multivariate random vector with mean 0 for null
#'  and mean \emph{mu} for alternatives. The second sample of \emph{n} \emph{m}-dimendional observations is drawn from with mean zero.
#'  The correlation matrix \emph{R} is block diagonal with a given \emph{block size} and with off-diagonal entries of \emph{correlation}.
#'
#'
#' @param m Number of tests
#' @param n Number of observations in each of the two samples
#' @param mu The mean under the alternative
#' @param prob The probability of alternatives
#' @param block_size The size of the blocks in the block-diagonal correlation matrix R.
#' @param rho The off-diagonal entries of th block matrices
#'
#' @return A list, with P-values (P), the pooled variance as the informative covariate (X), and the correlation matrix (Sigma).
#'
#' @importFrom Matrix Matrix bdiag
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixTests col_t_equalvar
#' @export
toy_example <- function(m=10000,n=15,mu=1,prob=0.05, block_size=100,rho=0.8){


  mat <- Matrix::Matrix((1-rho)*diag(block_size)+matrix(rho,block_size,block_size))
  lmat <- list()
  Z <- NULL
  for(i in 1:ceiling(m/block_size)){
    lmat <-c(lmat,mat)
    Z <- cbind(Z,mvtnorm::rmvnorm(n = 2*n, sigma=as.matrix(mat)))
  }

  Sigma <- Matrix::bdiag(lmat)
  Sigma <- Sigma[1:m,1:m]
  Z0 <- Z[,1:m]


  H <- rbinom(n = m,size = 1,prob = prob)
  H <- matrix(H,nrow=2*n,ncol=m,byrow = T)
  H[(n+1):(2*n)
    ,] <- 0
    Z <- data.frame(Z0+H*mu)
    X <- apply(Z,2,var)
    P <- matrixTests::col_t_equalvar(Z[1:n,],Z[(n+1):(2*n),])$pvalue
    return(list(P=P,X=X,Sigma=Sigma))

}
