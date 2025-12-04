#' @export
build_paths <- function(x, y, family, K, eps, delta, L = NULL) {
  
  #convert x to data frame and get variable names
  x <- as.data.frame(x)
  var_names <- colnames(x)
  if (is.null(var_names)) {
    var_names <- paste0("X", 1:ncol(x))
    colnames(x) <- var_names
  }
  
  #initialize
  frontiers <- list()
  aic_by_model <- list()
  
  #start with empty model S_0
  S_0 <- fit_empty_model(y, family)
  frontiers[[1]] <- data.frame(
    step = 0,
    variables = I(list(character(0))),
    aic = S_0$aic,
    model_id = model_id(character(0)),
    stringsAsFactors = FALSE
  )
  aic_by_model[[model_id(character(0))]] <- S_0$aic
  parent_best_aic <- S_0$aic
  
  #repeat for K steps
  for (k in 1:K) {
    parent_models <- frontiers[[k]]
    candidate_children <- list()
    
    for (i in 1:nrow(parent_models)) {
      parent_vars <- parent_models$variables[[i]]
      parent_aic <- parent_models$aic[i]
      
      available_vars <- setdiff(var_names, parent_vars)
      if (length(available_vars) == 0) next
      
      child_aics <- numeric(length(available_vars))
      
      for (j in seq_along(available_vars)) {
        child_vars <- c(parent_vars, available_vars[j])
        child_id <- model_id(child_vars)
        
        if (child_id %in% names(aic_by_model)) {
          child_aics[j] <- aic_by_model[[child_id]]
        } else {
          result <- fit_model(child_vars, x, y, family)
          child_aics[j] <- result$aic
          aic_by_model[[child_id]] <- result$aic
        }
      }
      
      best_child_aic <- min(child_aics)
      
      for (j in seq_along(available_vars)) {
        aic_improvement <- parent_aic - child_aics[j]
        within_delta <- (child_aics[j] - best_child_aic) <= delta
        improves_enough <- aic_improvement >= eps
        
        if (within_delta && improves_enough) {
          child_vars <- c(parent_vars, available_vars[j])
          candidate_children[[length(candidate_children) + 1]] <- list(
            variables = child_vars,
            aic = child_aics[j],
            model_id = model_id(child_vars),
            parent_id = parent_models$model_id[i]
          )
        }
      }
    }
    
    if (length(candidate_children) == 0) break
    
    children_df <- data.frame(
      step = k,
      variables = I(lapply(candidate_children, function(x) x$variables)),
      aic = sapply(candidate_children, function(x) x$aic),
      model_id = sapply(candidate_children, function(x) x$model_id),
      parent_id = sapply(candidate_children, function(x) x$parent_id),
      stringsAsFactors = FALSE
    )
    
    children_df <- children_df[!duplicated(children_df$model_id), ]
    
    current_best_aic <- min(children_df$aic)
    if ((parent_best_aic - current_best_aic) < eps) break
    
    if (!is.null(L) && nrow(children_df) > L) {
      children_df <- children_df[order(children_df$aic)[1:L], ]
    }
    
    frontiers[[k + 1]] <- children_df
    parent_best_aic <- current_best_aic
  }
  
  path_forest <- list(
    frontiers = frontiers,
    aic_by_model = aic_by_model,
    meta = list(
      family = family,
      K = K,
      eps = eps,
      delta = delta,
      L = L,
      n_steps = length(frontiers) - 1,
      n_vars = length(var_names),
      var_names = var_names,
      total_models_explored = length(aic_by_model)
    )
  )
  
  return(path_forest)
}
