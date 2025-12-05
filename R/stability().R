#' Compute Variable Stability Scores via Resampling
#'
#' Estimates the stability of variable selection by repeatedly building path forests
#' on resampled data and tracking how often each variable appears across models.
#' Higher stability scores indicate variables that are consistently selected across
#' different data samples, suggesting more robust importance.
#'
#' @param x A matrix or data frame of predictor variables. Column names are used
#'   as variable identifiers.
#' @param y A vector of response values. Should be numeric for linear regression
#'   or binary (0/1) for logistic regression.
#' @param B An integer specifying the number of resampling iterations. More iterations
#'   give more stable estimates but increase computation time. Default is 100.
#' @param resample A character string specifying the resampling method. Either
#'   "bootstrap" for sampling with replacement (n observations), or "subsample"
#'   for m-out-of-n subsampling without replacement. Default is "bootstrap".
#' @param m An integer specifying the subsample size when \code{resample = "subsample"}.
#'   Typically set to around 0.5*n to 0.8*n. Ignored if \code{resample = "bootstrap"}.
#' @param family A character string specifying the model family. Either "gaussian"
#'   for linear regression or "binomial" for logistic regression. Default is "gaussian".
#' @param K An integer specifying the maximum number of steps for path building
#'   in each resampling iteration. Default is 5.
#' @param eps A numeric threshold for minimum AIC improvement in path building.
#'   Default is 1e-6.
#' @param delta A numeric threshold for model competitiveness in path building.
#'   Default is 2.
#' @param L An integer specifying the maximum frontier size in path building.
#'   Controls computational cost in each iteration.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{pi}{A named numeric vector of stability scores for each variable,
#'       ranging from 0 to 1. Higher values indicate more stable selection.}
#'     \item{resample_proportions}{A B × p matrix where entry (b, j) is the
#'       proportion of models in resample b that included variable j. Useful
#'       for examining stability variability across resamples.}
#'   }
#'
#' @details
#' The stability estimation procedure works as follows:
#' \enumerate{
#'   \item For each of B resampling iterations:
#'     \itemize{
#'       \item Draw a resample of the data (bootstrap or subsample)
#'       \item Build a path forest using \code{\link{build_paths}}
#'       \item Extract all unique models across all frontiers
#'       \item For each variable, compute: (number of models containing it) / (total models)
#'     }
#'   \item Average these proportions across all B iterations to get stability score pi
#' }
#'
#' Stability scores near 1 indicate variables consistently selected across resamples,
#' while scores near 0 indicate rarely selected variables. The \code{resample_proportions}
#' matrix allows examination of variability in selection across resamples.
#'
#' Bootstrap tends to give higher variance estimates, while subsampling (with m < n)
#' can provide better finite-sample properties. Common choices are m ≈ 0.5n or 0.632n.
#'
#' @examples
#' # Linear regression example
#' set.seed(123)
#' n <- 100
#' x <- matrix(rnorm(n * 8), ncol = 8)
#' colnames(x) <- paste0("X", 1:8)
#' # True model uses X1, X2, X3
#' y <- x[,1] + 2*x[,2] - x[,3] + rnorm(n)
#'
#' # Compute stability with bootstrap
#' stab_boot <- stability(x, y, B = 50, resample = "bootstrap",
#'                        family = "gaussian", K = 5, eps = 0.5, 
#'                        delta = 2, L = 10)
#' print(sort(stab_boot$pi, decreasing = TRUE))
#'
#' # Compute stability with subsampling
#' stab_sub <- stability(x, y, B = 50, resample = "subsample", m = 50,
#'                       family = "gaussian", K = 5, eps = 0.5,
#'                       delta = 2, L = 10)
#' print(sort(stab_sub$pi, decreasing = TRUE))
#'
#' # Visualize stability scores
#' barplot(sort(stab_boot$pi, decreasing = TRUE),
#'         main = "Variable Stability Scores",
#'         ylab = "Stability (pi)",
#'         las = 2)
#' abline(h = 0.5, col = "red", lty = 2)  # Common threshold
#'
#' # Examine variability across resamples
#' boxplot(stab_boot$resample_proportions,
#'         main = "Selection Proportions Across Resamples",
#'         ylab = "Proportion",
#'         las = 2)
#'
#' # Use stability scores for model selection
#' paths <- build_paths(x, y, family = "gaussian", K = 5, 
#'                      eps = 0.5, delta = 2, L = 10)
#' plausible <- plausible_models(paths, stab_boot$pi, 
#'                               Delta = 2, tau = 0.6)
#'
#' @export
stability <- function(x, y, B = 100, resample = "bootstrap", m, family = "gaussian", K = 5, eps = 1e-6, delta = 2, L) {
  
  # remove 'AsIs' class if present
  if ("AsIs" %in% class(x)) {
    x <- unclass(x)
  }
  
  # initialize
  x <- as.data.frame(x)
  predictors <- colnames(x)
  p <- length(predictors)
  n <- length(y)
  Z <- matrix(0, nrow = B, ncol = p)
  colnames(Z) <- predictors
  
  # run bootstrap or subsampling
  for (b in seq_len(B)) {
    if (resample == "bootstrap") {
      idx <- sample(n, n, replace = TRUE)
    } else {  # m-out-of-n subsampling
      idx <- sample(n, m, replace = FALSE)
    }
    
    # build the path forest on the resampled data
    paths <- build_paths(
      x = x[idx, , drop = FALSE],
      y = y[idx],
      family = family,
      K = K,
      eps = eps,
      delta = delta,
      L = L
    )
    
    # collect unique model IDs across all frontiers
    all_model_ids <- character(0)
    all_models_list <- list()
    
    for (frontier in paths$frontiers) {
      if (nrow(frontier) > 0) {
        for (i in 1:nrow(frontier)) {
          mid <- frontier$model_id[i]
          if (!(mid %in% all_model_ids)) {
            all_model_ids <- c(all_model_ids, mid)
            all_models_list[[mid]] <- frontier$variables[[i]]
          }
        }
      }
    }
    
    n_models <- length(all_models_list)
    
    # update selection proportions
    if (n_models > 0) {
      feature_counts <- integer(p)
      names(feature_counts) <- predictors
      
      for (model_vars in all_models_list) {
        if (length(model_vars) > 0) {
          feature_counts[model_vars] <- feature_counts[model_vars] + 1
        }
      }
      
      Z[b, ] <- feature_counts / n_models
    }
  }
  
  # stability scores
  pi <- colMeans(Z)
  
  list(pi = pi, resample_proportions = Z)
}
