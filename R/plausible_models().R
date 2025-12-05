#' Select Plausible Models from Path Forest
#'
#' Filters models from a path forest based on AIC competitiveness, variable
#' stability, and similarity to identify a set of plausible final models.
#' This function helps narrow down the many models explored during path building
#' to a manageable set of high-quality candidates.
#'
#' @param path_forest A path forest object returned by \code{\link{build_paths}},
#'   containing frontiers of models and associated metadata.
#' @param pi A named numeric vector of stability scores for each variable, typically
#'   computed by \code{\link{compute_stability}}. Names should match variable names
#'   in the path forest.
#' @param Delta A numeric threshold for AIC competitiveness. Only models with
#'   AIC <= min(AIC) + Delta are retained. Larger values include more models.
#'   Default is 2.
#' @param tau A numeric threshold for average variable stability. Only models with
#'   mean stability score >= tau are retained. Values should be between 0 and 1.
#'   Default is 0.5.
#' @param threshold A numeric threshold for Jaccard similarity (0 to 1) used to
#'   remove near-duplicate models. If two models have Jaccard similarity >= threshold,
#'   the one with higher AIC is removed. If NULL, no similarity filtering is applied.
#'   Default is NULL.
#'
#' @return A data frame of plausible models, sorted by AIC (best first), with columns:
#'   \describe{
#'     \item{model_id}{Unique identifier for the model based on its variable set.}
#'     \item{variables}{List column containing character vectors of variable names
#'       in each model.}
#'     \item{aic}{AIC value for the model.}
#'     \item{pi_bar}{Average stability score across all variables in the model.}
#'     \item{n_vars}{Number of variables in the model.}
#'   }
#'
#' @details
#' The function applies a sequence of filters to identify plausible models:
#' \enumerate{
#'   \item \strong{Extract unique models}: Collects all unique models across all
#'     frontiers in the path forest.
#'   \item \strong{AIC filter}: Keeps only models within Delta AIC units of the
#'     best model.
#'   \item \strong{Stability filter}: Computes average stability (pi_bar) for each
#'     model and keeps only those with pi_bar >= tau.
#'   \item \strong{Similarity filter} (if threshold specified): Removes near-duplicate
#'     models using Jaccard similarity. For similar model pairs, keeps the one with
#'     lower AIC.
#' }
#'
#' Jaccard similarity between two models is calculated as:
#' \deqn{J(A,B) = \frac{|A \cap B|}{|A \cup B|}}
#' where A and B are the variable sets of two models.
#'
#' @examples
#' # Build paths
#' set.seed(123)
#' n <- 100
#' x <- matrix(rnorm(n * 5), ncol = 5)
#' colnames(x) <- paste0("Var", 1:5)
#' y <- x[,1] + 2*x[,2] + rnorm(n)
#'
#' paths <- build_paths(x, y, family = "gaussian", K = 5, eps = 0.5, delta = 2)
#'
#' # Compute stability (example - normally use compute_stability)
#' pi <- setNames(runif(5, 0.3, 0.9), colnames(x))
#'
#' # Select plausible models with default thresholds
#' plausible <- plausible_models(paths, pi)
#' print(plausible)
#'
#' # More stringent AIC and stability criteria
#' plausible_strict <- plausible_models(paths, pi, Delta = 1, tau = 0.7)
#'
#' # Remove similar models
#' plausible_diverse <- plausible_models(paths, pi, Delta = 2, tau = 0.5, 
#'                                       threshold = 0.7)
#'
#' # Examine the selected models
#' for (i in 1:nrow(plausible)) {
#'   cat("Model", i, ":", paste(plausible$variables[[i]], collapse = ", "), "\n")
#'   cat("  AIC:", round(plausible$aic[i], 2), 
#'       "  Avg Stability:", round(plausible$pi_bar[i], 2), "\n\n")
#' }
#'
#' @seealso 
#'   \code{\link{build_paths}} for creating the path forest,
#'   \code{\link{compute_stability}} for calculating variable stability scores
#'
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
