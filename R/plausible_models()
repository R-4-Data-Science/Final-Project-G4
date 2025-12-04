#' @export
plausible_models <- function(path_forest, pi, Delta = 2, tau = 0.5, threshold = NULL) {
  
  # init
  all_models <- list()
  all_model_ids <- character(0)
  
  # extract unique models
  for (frontier in path_forest$frontiers) {
    if (nrow(frontier) > 0) {
      for (i in seq_len(nrow(frontier))) {
        mid <- frontier$model_id[i]
        if (!(mid %in% all_model_ids)) {
          all_model_ids <- c(all_model_ids, mid)
          all_models[[mid]] <- list(
            model_id  = mid,
            variables = frontier$variables[[i]],
            aic       = frontier$aic[i]
          )
        }
      }
    }
  }
  
  # convert to df
  models_df <- data.frame(
    model_id  = sapply(all_models, `[[`, "model_id"),
    variables = I(lapply(all_models, `[[`, "variables")),
    aic       = sapply(all_models, `[[`, "aic"),
    stringsAsFactors = FALSE
  )
  
  # AIC filter
  aic_min <- min(models_df$aic)
  models_df <- models_df[models_df$aic <= aic_min + Delta, ]
  
  # average stability
  models_df$pi_bar <- sapply(models_df$variables, function(vars) {
    if (length(vars) == 0) return(0)
    mean(pi[vars])
  })
  
  # filter by stability threshold
  models_df <- models_df[models_df$pi_bar >= tau, ]
  
  # number of vars
  models_df$n_vars <- sapply(models_df$variables, length)
  
  # remove near-duplicates
  if (!is.null(threshold) && nrow(models_df) > 1) {
    keep_idx <- rep(TRUE, nrow(models_df))
    
    for (i in 1:(nrow(models_df)-1)) {
      if (!keep_idx[i]) next
      
      vars_i <- models_df$variables[[i]]
      
      for (j in (i+1):nrow(models_df)) {
        if (!keep_idx[j]) next
        
        vars_j <- models_df$variables[[j]]
        
        # Jaccard similarity
        intersection <- length(intersect(vars_i, vars_j))
        union <- length(union(vars_i, vars_j))
        jaccard <- if (union == 0) 0 else intersection / union
        
        # apply similarity threshold
        if (jaccard >= threshold) {
          if (models_df$aic[j] < models_df$aic[i]) {
            keep_idx[i] <- FALSE
            break
          } else {
            keep_idx[j] <- FALSE
          }
        }
      }
    }
    
    models_df <- models_df[keep_idx, ]
  }
  
  # sort by AIC
  models_df <- models_df[order(models_df$aic), ]
  rownames(models_df) <- NULL
  
  return(models_df)
}
