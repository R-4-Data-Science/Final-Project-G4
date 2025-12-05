#' Compute Confusion Matrix and Classification Metrics
#'
#' Calculates a confusion matrix and associated performance metrics for binary
#' classification models. Converts predicted probabilities to binary predictions
#' using a specified cutoff threshold, then computes standard diagnostic metrics.
#'
#' @param y_true A numeric vector of true binary outcomes (0 or 1).
#' @param y_prob A numeric vector of predicted probabilities, typically ranging
#'   from 0 to 1. Should be the same length as \code{y_true}.
#' @param cutoff A numeric threshold for converting probabilities to binary
#'   predictions. Probabilities >= cutoff are classified as 1, otherwise 0.
#'   Default is 0.5.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{confusion_matrix}{A 2x2 table showing the counts of True Negatives (TN),
#'       False Positives (FP), False Negatives (FN), and True Positives (TP).}
#'     \item{metrics}{A named numeric vector containing:
#'       \itemize{
#'         \item \strong{prevalence}: Proportion of positive cases in the data
#'         \item \strong{accuracy}: Overall proportion of correct predictions
#'         \item \strong{sensitivity}: True positive rate (recall); TP/(TP+FN)
#'         \item \strong{specificity}: True negative rate; TN/(TN+FP)
#'         \item \strong{precision}: Positive predictive value; TP/(TP+FP)
#'         \item \strong{FDR}: False discovery rate; FP/(TP+FP)
#'         \item \strong{DOR}: Diagnostic odds ratio; (TP/FN)/(FP/TN)
#'         \item \strong{N}: Total number of observations
#'         \item \strong{TP, TN, FP, FN}: Raw counts from confusion matrix
#'       }}
#'     \item{cutoff}{The cutoff threshold used for classification.}
#'   }
#'
#' @details
#' The function handles edge cases gracefully:
#' \itemize{
#'   \item If the confusion matrix is not fully populated (e.g., no predicted
#'     positives), it creates a complete 2x2 matrix with zeros where needed.
#'   \item Metrics involving division by zero return NA instead of causing errors.
#'   \item When both FN and FP are zero (perfect classifier), DOR returns Inf.
#' }
#'
#' The Diagnostic Odds Ratio (DOR) combines sensitivity and specificity into a
#' single metric. Higher values indicate better discrimination, with Inf representing
#' perfect classification.
#'
#' @examples
#' # Simulate some predictions
#' set.seed(123)
#' y_true <- rbinom(100, 1, 0.4)
#' y_prob <- plogis(rnorm(100, mean = ifelse(y_true == 1, 1, -1)))
#'
#' # Compute metrics with default cutoff
#' results <- confusion_metrics(y_true, y_prob)
#' print(results$confusion_matrix)
#' print(results$metrics)
#'
#' # Try different cutoff
#' results_strict <- confusion_metrics(y_true, y_prob, cutoff = 0.7)
#' print(results_strict$metrics["sensitivity"])
#'
#' # Compare multiple cutoffs
#' cutoffs <- seq(0.3, 0.7, by = 0.1)
#' accuracies <- sapply(cutoffs, function(c) {
#'   confusion_metrics(y_true, y_prob, cutoff = c)$metrics["accuracy"]
#' })
#' plot(cutoffs, accuracies, type = "b", 
#'      xlab = "Cutoff", ylab = "Accuracy")
#'
#' @seealso \code{\link{build_paths}} for building classification models,
#'   \code{\link{select_plausible_models}} for model selection
#'
#' @export
confusion_metrics <- function(y_true, y_prob, cutoff = 0.5) {
  #apply cutoff
  y_pred <- ifelse(y_prob >= cutoff, 1, 0)
  
  #create confusion matrix
  confusion <- table(Actual = y_true, Predicted = y_pred)
  
  #ensure 2x2 matrix
  if (nrow(confusion) < 2 || ncol(confusion) < 2) {
    #create matrix with zeros
    full_confusion <- matrix(0, nrow = 2, ncol = 2,
                             dimnames = list(Actual = c("0", "1"),
                                             Predicted = c("0", "1")))
    #fill in observed values
    for (i in rownames(confusion)) {
      for (j in colnames(confusion)) {
        full_confusion[i, j] <- confusion[i, j]
      }
    }
    confusion <- full_confusion
  }
  
  #extract counts
  TN <- confusion["0", "0"]
  FP <- confusion["0", "1"]
  FN <- confusion["1", "0"]
  TP <- confusion["1", "1"]
  
  N <- TN + FP + FN + TP
  
  #compute metrics
  prevalence <- (TP + FN) / N
  accuracy <- (TP + TN) / N
  
  #handle division by zero
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  FDR <- if ((TP + FP) > 0) FP / (TP + FP) else NA
  
  #DOR
  DOR <- if (FN > 0 && FP > 0 && TN > 0) {
    (TP / FN) / (FP / TN)
  } else if (FN == 0 && FP == 0) {
    Inf  # Perfect classifier
  } else {
    NA
  }
  
  #metrics vector
  metrics <- c(
    prevalence = prevalence,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    FDR = FDR,
    DOR = DOR,
    N = N,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  )
  
  #return list
  result <- list(
    confusion_matrix = confusion,
    metrics = metrics,
    cutoff = cutoff
  )
  return(result)
}
