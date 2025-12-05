\name{confusion_metrics}
\alias{confusion_metrics}
\title{Compute Confusion Matrix and Classification Metrics}
\usage{
confusion_metrics(y_true, y_prob, cutoff = 0.5)
}
\arguments{
\item{y_true}{A numeric vector of true binary outcomes (0 or 1).}

\item{y_prob}{A numeric vector of predicted probabilities, typically ranging
from 0 to 1. Should be the same length as \code{y_true}.}

\item{cutoff}{A numeric threshold for converting probabilities to binary
predictions. Probabilities >= cutoff are classified as 1, otherwise 0.
Default is 0.5.}
}
\value{
A list with three components:
\describe{
\item{confusion_matrix}{A 2x2 table showing the counts of True Negatives (TN),
False Positives (FP), False Negatives (FN), and True Positives (TP).}
\item{metrics}{A named numeric vector containing prevalence, accuracy,
sensitivity, specificity, precision, FDR, DOR, N, TP, TN, FP, and FN.}
\item{cutoff}{The cutoff threshold used for classification.}
}
}
\description{
Calculates a confusion matrix and associated performance metrics for binary
classification models. Converts predicted probabilities to binary predictions
using a specified cutoff threshold, then computes standard diagnostic metrics.
}
\details{
The function handles edge cases gracefully. If the confusion matrix is not 
fully populated, it creates a complete 2x2 matrix with zeros where needed.
Metrics involving division by zero return NA instead of causing errors.

The Diagnostic Odds Ratio (DOR) combines sensitivity and specificity into a
single metric. Higher values indicate better discrimination, with Inf representing
perfect classification.
}
\examples{
# Simulate some predictions
set.seed(123)
y_true <- rbinom(100, 1, 0.4)
y_prob <- plogis(rnorm(100, mean = ifelse(y_true == 1, 1, -1)))

# Compute metrics with default cutoff
results <- confusion_metrics(y_true, y_prob)
print(results$confusion_matrix)
print(results$metrics)

# Try different cutoff
results_strict <- confusion_metrics(y_true, y_prob, cutoff = 0.7)
print(results_strict$metrics["sensitivity"])
}
\seealso{
\code{\link{build_paths}} for building classification models,
\code{\link{plausible_models}} for model selection
}
