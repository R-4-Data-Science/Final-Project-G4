\name{stability}
\alias{stability}
\title{Compute Variable Stability Scores via Resampling}
\usage{
stability(x, y, B = 100, resample = "bootstrap", m, 
          family = "gaussian", K = 5, eps = 1e-06, delta = 2, L)
}
\arguments{
\item{x}{A matrix or data frame of predictor variables.}

\item{y}{A vector of response values.}

\item{B}{Number of resampling iterations. Default is 100.}

\item{resample}{Resampling method: "bootstrap" or "subsample". Default is "bootstrap".}

\item{m}{Subsample size when resample = "subsample".}

\item{family}{Model family: "gaussian" or "binomial". Default is "gaussian".}

\item{K}{Maximum number of steps for path building. Default is 5.}

\item{eps}{Minimum AIC improvement threshold. Default is 1e-6.}

\item{delta}{Model competitiveness threshold. Default is 2.}

\item{L}{Maximum frontier size in path building.}
}
\value{
A list with two components:
\describe{
\item{pi}{A named numeric vector of stability scores for each variable (0 to 1).}
\item{resample_proportions}{A B × p matrix of selection proportions across resamples.}
}
}
\description{
Estimates the stability of variable selection by repeatedly building path forests
on resampled data and tracking how often each variable appears across models.
}
\details{
For each of B resampling iterations:
\enumerate{
\item Draw a resample of the data (bootstrap or subsample)
\item Build a path forest
\item Extract all unique models
\item Compute proportion of models containing each variable
}

The stability score pi is the average of these proportions across all iterations.
}
\examples{
# Linear regression example
set.seed(123)
n <- 100
x <- matrix(rnorm(n * 8), ncol = 8)
colnames(x) <- paste0("X", 1:8)
y <- x[,1] + 2*x[,2] - x[,3] + rnorm(n)

# Compute stability with bootstrap
stab_boot <- stability(x, y, B = 50, resample = "bootstrap",
                      family = "gaussian", K = 5, eps = 0.5, 
                      delta = 2, L = 10)
print(sort(stab_boot$pi, decreasing = TRUE))
}
\seealso{
\code{\link{build_paths}} for the path building algorithm,
\code{\link{plausible_models}} for using stability scores
}
