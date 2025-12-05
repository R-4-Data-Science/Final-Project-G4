#' Build Model Selection Paths Using Multi-Path Algorithm
#'
#' Constructs a forest of model selection paths by iteratively adding variables
#' to models based on AIC improvement criteria. At each step, the algorithm explores
#' multiple promising paths simultaneously, keeping models that improve sufficiently
#' and are competitive with the best model at that step.
#'
#' @param x A matrix or data frame of predictor variables. Column names will be used
#'   as variable names; if missing, variables will be named X1, X2, etc.
#' @param y A vector of response values. Should be numeric for linear regression
#'   or binary (0/1) for logistic regression.
#' @param family A character string specifying the model family. Must be either
#'   "gaussian" for linear regression or "binomial" for logistic regression.
#' @param K An integer specifying the maximum number of steps (variables to add).
#'   The algorithm will terminate early if no models meet the improvement criteria.
#' @param eps A numeric threshold for minimum AIC improvement required to add a
#'   variable. Variables must improve parent model AIC by at least this amount.
#' @param delta A numeric threshold for competitiveness. Child models are kept if
#'   their AIC is within delta of the best child model at that step.
#' @param L An integer specifying the maximum number of models to retain at each
#'   step (frontier size). If NULL, all qualifying models are kept. Use this to
#'   control computational cost for high-dimensional problems.
#'
#' @return A list with class "path_forest" containing:
#'   \describe{
#'     \item{frontiers}{A list of data frames, one per step, containing the models
#'       at each frontier. Each data frame has columns: step, variables (list column),
#'       aic, model_id, and parent_id.}
#'     \item{aic_by_model}{A named list storing AIC values for all explored models,
#'       indexed by model_id to avoid redundant computation.}
#'     \item{meta}{A list of metadata including family, K, eps, delta, L, number of
#'       steps completed, number of variables, variable names, and total models explored.}
#'   }
#'
#' @details
#' The algorithm starts with an empty model (intercept only) and iteratively builds
#' a forest of model selection paths. At each step k:
#' \enumerate{
#'   \item For each model in the current frontier, evaluate adding each available variable
#'   \item Keep child models that: (a) improve parent AIC by at least eps, AND
#'     (b) are within delta AIC of the best child at this step
#'   \item Stop if no models meet criteria or K steps are reached
#' }
#'
#' The multi-path approach explores multiple promising directions simultaneously,
#' unlike traditional stepwise selection which follows a single path.
#'
#' @examples
#' # Linear regression example
#' set.seed(123)
#' n <- 100
#' x <- matrix(rnorm(n * 5), ncol = 5)
#' colnames(x) <- paste0("Var", 1:5)
#' y <- x[,1] + 2*x[,2] + rnorm(n)
#'
#' paths <- build_paths(x, y, family = "gaussian", K = 5, eps = 0.5, delta = 2)
#'
#' # Logistic regression example
#' y_binary <- rbinom(n, 1, plogis(x[,1] + x[,2]))
#' paths_logit <- build_paths(x, y_binary, family = "binomial", 
#'                            K = 5, eps = 0.5, delta = 2, L = 10)
#'
#' @seealso \code{\link{compute_stability}} for analyzing path stability,
#'   \code{\link{select_plausible_models}} for choosing final models
#'
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
