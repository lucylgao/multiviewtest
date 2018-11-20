#' Pseudo likelihood ratio test for independence of clusterings
#'
#'
#' Implements the optimization algorithm for solving equation 2.8)
#' in Section 2.3.2 of "Are Clusterings of Multiple Data Views Independent?"
#' 
#'
#' @param logphi1 log(phi1), where the n x K1 matrix phi1 is defined in equation 2.9
#' @param logphi2 log(phi1), where the n x K2 matrix phi2 is defined in equation 2.9
#' use for the permutation procedure. The default number is 200.
#' @param row K1-vector containing the estimated mixture component probabilities on view 1
#' @param col K2-vector containing the estimated mixture component probabilities on view 2
#' @param stepsize Fixed step size to be used in the optimization. See Supplement B of
#'  "Are Clusterings of Multiple Data Views Independent?" for details.
#' @param maxiter Maximum number of iterations to run. The default maximum is 1000.
#' @param Pi An optional initial value for the optimization algorithm. The default
#' value is matrix(1, K1, K2). 
#' @import matrixStats
#' @export
#' @return
#' A list containing the following output components:
#' \item{Pi.est}{The maximizer}
#' \item{obj}{The objective value.}
#' @examples
#' library(mclust)
#' set.seed(1)
#' n <- 50
#' sig <- 2
#' p <- 2
#' K <- 3
#' mu1 <- cbind(c(2, 0), c(0, 2),  c(2, -2), c(-2, 0), c(0, -2), c(-2, 2))
#' mu2 <- cbind(c(-2, 0), c(0, -2), c(-2, 2), c(2, 0), c(0, 2), c(2, -2))
#' # Generate two-view data where the clusters are identical, 
#' # so Pi = I/3.
#' cl <- sample(1:K, n, replace=TRUE)
#' x <- list(matrix(sig* rnorm(n*p), n, p) + t(mu1)[cl, ],
#' matrix(sig * rnorm(n*p), n, p) + t(mu2)[cl, ])
#'
#'# Clusters each view
#'EM.View1 <- Mclust(x[[1]], G=3, modelNames=c("EII"))
#'EM.View2 <- Mclust(x[[2]], G=3, modelNames=c("EII"))
#'
#'# Constructs phi1 and phi2 matrices
#'logphi1 <-  cdens("EII", x[[1]], logarithm=TRUE, EM.View1$parameters)
#'logphi2 <-  cdens("EII", x[[2]], logarithm=TRUE, EM.View2$parameters)
#' 
#' pi1.est <- EM.View1$parameters$pro
#' pi2.est <- EM.View2$parameters$pro 
#' 
#' # Solves equation 2.8
#' Pi.est <- optimize_over_pi_clust(logphi1, logphi2,
#'  pi1.est, pi2.est, step=0.001)
#'  Pi.est$Pi
#' 
#'
#' @references
#' Gao, L.L., Bien, J., Witten, D. (2018) Are Clusterings of Multiple Data Views Independent?
#' submitted to Biostatistics.
#'

optimize_over_pi_clust <- function(logphi1, logphi2, row,col,
                             maxiter = 1000, stepsize = 0.01, Pi = NULL) {
  stopifnot(nrow(logphi1) == nrow(logphi2))
  n <- nrow(logphi1)
  K1 <- ncol(logphi1)
  K2 <- ncol(logphi2)
  if (is.null(Pi)) Pi <- matrix(1 / (K1 * K2), K1, K2)
  obj <- rep(NA, maxiter)

  # Creates a 3-D array of dimension n x K1 x K2
  # where the element in position ikk' contains
  # the sum-normalized phi1[i, k]*the sum-normalized phi2[i, k']

  # Uses the exp normalize trick for stability:
  # phi1/rowSums(phi1) = exp(logphi1)/rowSums(exp(logphi1))
  # = exp(logphi1 - rowMaxes(exp(logphi1)))/rowSums(exp(logphi1 - rowMaxes(exp(logphi1))))
  #normalized_phi1 <- phi1/rowSums(phi1)
  #normalized_phi2 <- phi2/rowSums(phi2)
  logphi1_standard <- logphi1 - rowMaxs(logphi1)
  normalized_phi1 <- exp(logphi1_standard)/rowSums(exp(logphi1_standard))
  logphi2_standard <- logphi2 - rowMaxs(logphi2)
  normalized_phi2 <- exp(logphi2_standard)/rowSums(exp(logphi2_standard))
  phi_3D <- array(rep(normalized_phi1, K2), c(n, K1, K2))*
    aperm(array(rep(t(normalized_phi2),each=K1), c(K1, K2, n)), c(3, 1, 2))

  for (j in seq(maxiter)) {
    # Creates a 3-D array of dimension n x K1 x K2
    # where the element in position ikk' contains
    # the sum-normalized phi1[i, k]*the sum-normalized phi2[i, k']*Pi[k, k']
    phipi_3D <- array(rep(Pi, each=n), dim=c(n, K1, K2))*phi_3D

    # The components of the interior sum in the objective function
    objVect <- rowSums(phipi_3D)
    obj[j] <- -sum(log(objVect))
    
    # If the relative tolerance is less than 1e-10, terminate
    if(abs(obj[j] - obj[j-1])/(min(obj[j], obj[j-1])) < 1e-10 && j > 1) {
      break
    }

    # Computes G, K1 x K2 matrix
    G <- crossprod(normalized_phi1, normalized_phi2/objVect)
    
    # Computes M, K1 x K2 matrix. We use Sinkhorn-Knopp
    # To scale the rows and columns to Pi_1 hat
    # and Pi_2 hat respectively.
    M <- Pi * exp(stepsize * G - 1)

    if(!is.double(sum(M))) stop("Gradient calculation failed")    
    
    # Sinkhorn-Knopp algorithm from
    # Sinkhorn Distances: Lightspeed Computation of Optimal Transport
    # Cuturi (2013)
    v = rep(1, K2)
    u = rep(1, K1)

    K = M/row
    for(i in 1:100) {
      uprev=u
      u = 1/K%*%(col/t(M)%*%u)

      # Convergence criteria
      if(abs(uprev - u)/min(uprev, u) < 1e-6 && j > 1) {
        break
      }
    }

    v = as.numeric(col/(t(M)%*%u))

    # Update Pi
    Pi <- diag(as.numeric(u))%*%M%*%diag(v)
  }

  obj <- obj - (sum(log(rowSums(exp(logphi1_standard)))) +
                  sum(log(rowSums(exp(logphi2_standard)))) + sum(rowMaxs(logphi1)) +
                  sum(rowMaxs(logphi2)))

  list(Pi = Pi, obj = obj)

}
