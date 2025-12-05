\name{plausible_models}
\alias{plausible_models}
\title{Select Plausible Models from Path Forest}
\usage{
plausible_models(path_forest, pi, Delta = 2, tau = 0.5, threshold = NULL)
}
\arguments{
\item{path_forest}{A path forest object returned by \code{\link{build_paths}},
containing frontiers of models and associated metadata.}

\item{pi}{A named numeric vector of stability scores for each variable, typically
computed by \code{\link{stability}}. Names should match variable names
in the path forest.}

\item{Delta}{A numeric threshold for AIC competitiveness. Only models with
AIC <= min(AIC) + Delta are retained. Default is 2.}

\item{tau}{A numeric threshold for average variable stability. Only models with
mean stability score >= tau are retained. Values should be between 0 and 1.
Default is 0.5.}

\item{threshold}{A numeric threshold for Jaccard similarity (0 to 1) used to
remove near-duplicate models. If NULL, no similarity filtering is applied.
Default is NULL.}
}
\value{
A data frame of plausible models, sorted by AIC, with columns:
model_id, variables (list column), aic, pi_bar, and n_vars.
}
\description{
Filters models from a path forest based on AIC competitiveness, variable
stability, and similarity to identify a set of plausible final models.
}
\details{
The function applies a sequence of filters:
\enumerate{
\item Extract unique models from all frontiers
\item AIC filter: Keep only models within Delta AIC units of the best model
\item Stability filter: Keep only models with average stability >= tau
\item Similarity filter: Remove near-duplicates using Jaccard similarity
}

Jaccard similarity between two models is |A ∩ B| / |A ∪ B|.
}
\examples{
# Build paths
set.seed(123)
n <- 100
x <- matrix(rnorm(n * 5), ncol = 5)
colnames(x) <- paste0("Var", 1:5)
y <- x[,1] + 2*x[,2] + rnorm(n)

paths <- build_paths(x, y, family = "gaussian", K = 5, eps = 0.5, delta = 2)

# Compute stability (example)
pi <- setNames(runif(5, 0.3, 0.9), colnames(x))

# Select plausible models
plausible <- plausible_models(paths, pi)
print(plausible)
}
\seealso{
\code{\link{build_paths}} for creating the path forest,
\code{\link{stability}} for calculating stability scores
}
