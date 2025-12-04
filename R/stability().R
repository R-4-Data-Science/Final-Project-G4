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
