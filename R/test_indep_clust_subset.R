#' Pseudo likelihood ratio test for independence of clusterings
#'
#'
#' Implements the pseudo likelihood ratio test described in Section 3 of
#' "Are Clusterings of Multiple Data Views Independent?"
#' for testing whether clusterings of two data views are independent, 
#' when observations can be  removed in one view but not the other. 
#' Allows for fitting model-based clusterings on the full data views, 
#' and comparing these clusterings on subsets of the data views where 
#' observations are not removed in both views. 
#'
#' @param x Multi-view data with two views; a list of two numeric vectors
#'  (in the case of univariate data) or matrices containing the two data views.
#'  In matrix format, rows correspond to observations and columns correspond to variables.
#' @param B An integer specifying the number of permutations to
#' use for the permutation procedure. The default number is 200.
#' @param model1 A character string indicating the model
#' to be fitted for Gaussian model-based clustering in view 1 using
#' the function \code{\link[mclust]{Mclust}}. The default is \code{"EII"}
#' (spherical, equal volume).
#' The help file for \code{\link[mclust]{mclustModelNames}} describes
#' the available model options.
#' @param model2 A character string indicating indicating the model
#' to be fitted for Gaussian model-based clustering in view 1 using
#' the function \code{\link[mclust]{Mclust}}. The default is \code{"EII"}
#' (spherical, equal volume). The help file for
#'  \code{\link[mclust]{mclustModelNames}} describes the available model options.
#' @param subset1 A numeric vector indicating the rows of the matrix containing 
#' the first data view that correspond to observations which have not been 
#' removed in the second data view.
#' @param subset2 A numeric vector indicating the rows of the matrix containing 
#' the second data view that correspond to observations which have not been removed 
#' in the first data view.
#' @param step A numeric value containing the fixed step size to be used in the optimization
#' algorithm for estimating Pi. The default step size is 0.001. See Appendix B 
#' in the Supplementary Materials of "Are Clusterings of Multiple Data Views Independent?" 
#' for details.
#' @param maxiter A numeric value containing the maximum number of iterations to run in
#'  the optimization algorithm. The default maximum is 1000.
#' @param init1 An optional argument containing the model to be fitted in the
#' hierarchical clustering initialization in Gaussian model-based clustering
#' in view 1. The default is \code{"VVV"} (ellipsoidal, varying volume,
#'  shape, and orientation).
#' The help file for \code{\link[mclust]{hc}} describes
#' the available model options.
#' @param init2 An optional argument containing the model to be fitted in the
#' hierarchical clustering initialization in Gaussian model-based clustering
#' in view 2. The default is \code{"VVV"} (ellipsoidal, varying volume,
#'  shape, and orientation). 
#' The help file for \code{\link[mclust]{hc}} describes
#' the available model options.
#' @param K1 An optional argument containing the number of clusters in View 1.
#'  If left out, then the number of clusters is chosen with BIC as described in
#'  Section 2.3.3 of "Are Clusterings of Multiple Data Views Independent?"
#' @param K2 An optional argument containing the number of clusters in View 2.
#'  If left out, then the number of clusters is chosen with BIC as described in
#'  Section 2.3.3 of "Are Clusterings of Multiple Data Views Independent?"
#' @import mclust
#' @export
#' @return
#' A list containing the following output components:
#' \item{K1}{The number of clusters in view 1}
#' \item{K2}{The number of clusters in view 2}
#' \item{Pi.est}{The estimated Pi matrix}
#' \item{PLRstat}{The pseudo likelihood ratio test statistic}
#' \item{pval}{The p-value}
#' \item{modelfit1}{The object of class '\code{Mclust}' corresponding to the model-based clustering fitted
#' in View 1; contains eg. estimated parameters and cluster assignments. The help file for \code{\link[mclust]{Mclust}} describes the components of the
#' object.}
#' \item{modelfit2}{The object of class '\code{Mclust}' corresponding to the model-based clustering fitted
#' in View 2; contains eg. estimated parameters and cluster assignments.
#' The help file for \code{\link[mclust]{Mclust}} describes the components of the
#' object.}
#' @examples
#' set.seed(1)
#' n <- 50
#' sig <- 2
#' p <- 2
#' K <- 3
#' mu1 <- cbind(c(2, 0), c(0, 2),  c(2, -2), c(-2, 0), c(0, -2), c(-2, 2))
#' mu2 <- cbind(c(-2, 0), c(0, -2), c(-2, 2), c(2, 0), c(0, 2), c(2, -2))
#' # Generates two-view data where the clusters are independent.
#' x1 <- list(matrix(sig* rnorm(n*p), n, p) + t(mu1)[sample(1:K, n, replace=TRUE), ],
#'         matrix(sig * rnorm(n*p), n, p) + t(mu2)[sample(1:K, n, replace=TRUE), ])
#' # Generate two-view data where the clusters are identical.
#' n <- 71
#' cl <- sample(1:K, n, replace=TRUE)
#' x2 <- list(matrix(sig* rnorm(n*p), n, p) + t(mu1)[cl, ],
#' matrix(sig * rnorm(n*p), n, p) + t(mu2)[cl, ])
#'
#' # Run the function on independent data views; we do not reject the null hypothesis.
#' # Note: Will take a few seconds to run
#'  # By default, not specifying K1 and K2 means the number of clusters
#'  # to use in the test in each view is chosen via BIC.
#'  # Covariance matrix model specified is shared sigma^2 I covariance matrix in view 1
#'  # and shared diagonal covariance matrix in view 2.
#'  # B specifies the number of permutations to do for the permutation test.
#'  # Covariance matrix model specified for initialization
#'  # is shared sigma^2 I covariance matrix in view 1
#'  # Estimates Gaussian mixture model parameters on x1[[1]] and x1[[2]], 
#'  # and compares the estimated clusterings on the subsetted data 
#'  # x1[[1]][2:48, ] and x1[[1]][2:48, ].
#' indep1 <- test_indep_clust_subset(x1,model1="EII", model2="EEI",
#' subset1=2:48, subset2=2:48, init1="EII", B=52)
#' # The estimated cluster parameters in view 1
#' indep1$modelfit1$parameters
#' # The cluster assignments in view 2
#' indep1$modelfit2$classification
#' 
#' # Run the function on identical data views; we reject the null hypothesis
#' # Note: Will take a few seconds to run
#' # Additionally specify the number of clusters in each view to use in the test
#' # Covariance matrix model specified is shared covariance matrix in view 1
#' # and shared diagonal covariance matrix in view 2.
#' # See mclust documentation for more covariance model specification options.
#' identical2 <- test_indep_clust_subset(x2,model1="EEE", model2="EEI", 
#'                                subset1=1:70, subset2=1:70, 
#'                                K1=2, K2=3, B=51)
#' # P-value
#' identical2$pval
#'
#' @references
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016)
#' mclust 5: clustering, classification and density estimation
#' using Gaussian finite mixture models, The R Journal, 8/1, pp. 205-233.
#'
#' Fraley C. and Raftery A. E. (2002) Model-based clustering,
#' discriminant analysis and density estimation,
#' Journal of the American Statistical Association, 97/458, pp. 611-631.
#'
#' Gao, L.L., Bien, J., Witten, D. (2018) Are Clusterings of Multiple Data Views Independent?
#' submitted to Biostatistics.
#'
test_indep_clust_subset <- function(x,  model1="EII", model2="EII",
                              K1=NULL, K2=NULL, 
                              subset1, subset2,
                             init1=NULL, init2=NULL, 
                             B=200, step=0.001, maxiter=1000) {

  if(typeof(x) != "list") {
    stop("x is not in the right format; should be a list of two matrices
         (or vectors, for univariate data)")
  } else {
    if(length(x) != 2) {
      stop("x is not in the right format; should be a list of two numeric matrices
         (or vectors, for univariate data)")
    }  else {
      if(!is.numeric(x[[1]]) || !is.numeric(x[[2]])) {
        stop("x is not in the right format; should be a list of two numeric matrices
         (or vectors, for univariate data)")
      } else {
        if(is.matrix(x[[1]]) != is.matrix(x[[2]])) {
          stop("x is not in the right format; should be a list of two numeric matrices
         (or vectors, for univariate data)")
        } 
      }
    }
  }
  
  if(any(is.na(x[[1]])) || any(is.na(x[[2]])))
    stop(paste("views cannot contain any missing values"))
  
  if(!is.numeric(subset1) || !is.numeric(subset2)) {
    stop("subset is not in the right format; should be a numeric vector")    
  } else { 
    if(sum(c(subset1 <= 0, subset2 <= 0)) > 0) { 
      stop("cannot provide a subset with a negative index")  
    } else { 
      if(length(subset1) != length(subset2)){
        stop("the two subsets should contain the same number of observations")
      }
    }
  }
  if(!model1 %in% c("E", "V", mclust.options("emModelNames")))
    stop(paste("Did not pass a valid model for Mclust in View 1") )
  
  if(model1 %in% c("E", "V") && is.matrix(x[[1]]))
    stop(paste(model1, "is a univariate model for Mclust, cannot be used for
               multivariate data") )
  
  if(model1 %in% mclust.options("emModelNames") && !is.matrix(x[[1]]))
    stop(paste(model1, "is a multivariate model for Mclust, cannot be used for
               univariate data"))
  
  if(!model2 %in% c("E", "V", mclust.options("emModelNames")))
    stop(paste("Did not pass a valid model for Mclust in View 2") )
  
  if(model2 %in% c("E", "V") && is.matrix(x[[2]]))
    stop(paste(model2, "is a univariate model for Mclust, cannot be used for
               multivariate data"))
  if(model2 %in% mclust.options("emModelNames") && !is.matrix(x[[2]]))
    stop(paste(model2, "is a multivariate model for Mclust, cannot be used for
               univariate data") )
  
  
  if(is.null(init1)) { 
    init1 <- "VVV"  
  } else { 
    if(!init1 %in% c("E", "V", mclust.options("hcModelNames")))
      stop(paste("Did not pass a valid initialization model for Mclust in View 1") )
    
    if(init1 %in% c("E", "V") && is.matrix(x[[1]]))
      stop(paste(init1, "is a univariate initialization model for Mclust, 
                 cannot be used for multivariate data") )
    
    if(init1 %in% mclust.options("hcModelNames") && !is.matrix(x[[1]]))
      stop(paste(init1, "is a multivariate initialization model for Mclust, 
                 cannot be used for univariate data"))
    
    agree1 <- ((init1 %in% c("E", "V")) == (model1 %in% c("E", "V"))) 
    if(!agree1) { 
      stop("In view 1, initialization model is for univariate data, 
           but model is for multivariate data, or 
           initialization model is for multivariate data,
           but model is for univariate data")
    }    
    }
  
  if(is.null(init2)) { 
    init2 <-  "VVV"  
  } else { 
    if(!init2 %in% c("E", "V", mclust.options("hcModelNames")))
      stop(paste("Did not pass a valid initialization model for Mclust in View 2") )
    
    if(init2 %in% c("E", "V") && is.matrix(x[[2]]))
      stop(paste(init2, "is a univariate initialization model for Mclust, 
                 cannot be used for multivariate data"))
    if(init2 %in% mclust.options("hcModelNames") && !is.matrix(x[[2]]))
      stop(paste(init2, "is a univariate initialization model for Mclust,
                 cannot be used for univariate data") )
    agree2 <- ((init2 %in% c("E", "V")) == (model2 %in% c("E", "V"))) 
    if(!agree2) { 
      stop("In view 2, initialization model is for univariate data, 
           but model is for multivariate data, or 
           initialization model is for multivariate data,
           but model is for univariate data")
    }
 }
  
  if(step <=0 ) stop('Step size for optimization algorithm should be > 0')
  if(B < 1) stop('Number of permutations should be > 0')
  if(B <= 50) warning('Number of permutation iterations is specified to be small;
                     p-value approximation will likely be imprecise')
  if(maxiter %% 1 != 0|| maxiter <= 0) stop('Maximum number of iterations should be
                                                a positive integer')
  
  # Model-based clustering of each view
  if(is.null(K1)) {
    EM.View1 <- mclust::Mclust(x[[1]], G=2:9, modelNames=c(model1),
                               initialization=list(hcPairs=hc(x[[1]],
                                                              modelName=init1)))
    if(is.null(EM.View1)) stop("The model specified for model-based clustering
                               could not be fitted in view 1")
    K1 <- EM.View1$G
  } else {
    if(!is.numeric(K1)) stop('Number of clusters must be a positive integer in View 1')
    if(K1 %% 1 != 0 || K1 < 0) stop('Number of clusters must be a positive integer in View 1')
    if(K1 < 1) stop('Number of clusters must be greater than 1 in View 1')
    
    
    EM.View1 <- mclust::Mclust(x[[1]], K1, modelNames=c(model1),
                               initialization=list(hcPairs=hc(x[[1]],
                                                              modelName=init1)))
    if(is.null(EM.View1)) stop("The model specified for model-based clustering
                               could not be fitted in view 1")
  }
  
  if(is.null(K2)) {
    EM.View2 <- mclust::Mclust(x[[2]], G=2:9, modelNames=c(model2),
                               initialization=list(hcPairs=hc(x[[2]],
                                                              modelName=init2)))
    if(is.null(EM.View2)) stop("The model specified for model-based clustering
                               could not be fitted in view 1")
    K2 <- EM.View2$G
  } else {
    if(!is.numeric(K2)) stop('Number of clusters must be a positive integer in View 2')
    if(K2 %% 1 != 0 || K2 < 0) stop('Number of clusters must be a positive integer in View 2')
    if(K2 < 1) stop('Number of clusters must be greater than 1 in View 2')
    EM.View2 <- mclust::Mclust(x[[2]], K2, modelNames=c(model2),
                               initialization=list(hcPairs=hc(x[[2]],
                                                              modelName=init2)))
    if(is.null(EM.View2)) stop("The model specified for model-based clustering
                               could not be fitted in view 2")
  }
  
  
  EM.View1.param <- EM.View1$parameters
  EM.View2.param <- EM.View2$parameters
  
  if(is.matrix(x[[1]])) { 
  # Density matrices for each mixture component
    logphi1 <-  mclust::cdens(model1, x[[1]][subset1, ], logarithm=TRUE, EM.View1.param)
  } else { 
    logphi1 <-  mclust::cdens(model1, x[[1]][subset1], logarithm=TRUE, EM.View1.param)
  }
  
  if(is.matrix(x[[1]])) { 
    logphi2 <- mclust::cdens(model2, x[[2]][subset2, ], logarithm=TRUE, EM.View2.param)
  } else { 
    logphi2 <- mclust::cdens(model2, x[[2]][subset2], logarithm=TRUE, EM.View2.param)
  }
    
  n <- nrow(logphi1)
  pi1.est <- EM.View1.param$pro
  pi2.est <- EM.View2.param$pro
  # Optimization to estimate Pi
  Pi.est <- optimize_over_pi_clust(logphi1, logphi2,
                             pi1.est, pi2.est,
                             stepsize=step, maxiter=maxiter)

  # Compute pseudo likelihood ratio test statistic
  nullLogLik <- EM.View1$loglik + EM.View2$loglik
  log.Lambda = -min(Pi.est$obj, na.rm=T) - nullLogLik
  PLRstat <- 2*log.Lambda

  cat("Computing p-value: ")
  # Permutation procedure
  PLRstat.perm <- rep(0, B)
  for(b in 1:B) {

    if(b %% 10^round(log10(B) - 1) == 0) {
      cat(round(b/B*100))
      cat("% of permutations done ... ")
    }

    if(b == 1001) {
      earlyPval <- mean(ifelse(PLRstat.perm[1:1000] >= PLRstat, 1, 0))
      if(earlyPval >= 0.1) {
        PLRstat.perm <- PLRstat.perm[1:1000]
        cat("stopped permuting at 1000 permutations (calculated permutation p-value >= 0.1)")
        cat("\n")
        cat("\n")
        break
      }
    }
    Pi.est.perm <- optimize_over_pi_clust(logphi1, logphi2[sample(1:n, replace=F), ],
                                    pi1.est, pi2.est, stepsize=step)
    log.Lambda.perm <- -min(Pi.est.perm$obj, na.rm=T) - nullLogLik
    PLRstat.perm[b] <- 2*log.Lambda.perm
  }
  cat("\n")
  cat("\n")

  pval.perm <- mean(ifelse(PLRstat.perm >= PLRstat, 1, 0))

  # Likelihood Ratio statistic
  return(list(PLRstat=PLRstat,
              Pi.est=Pi.est$Pi,
              K1=K1,
              K2=K2,
              pval=pval.perm,
              modelfit1=EM.View1,
              modelfit2=EM.View2))
}
