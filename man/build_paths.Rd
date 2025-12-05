\name{build_paths}
\alias{build_paths}
\title{Build Model Selection Paths Using Multi-Path Algorithm}
\usage{
build_paths(x, y, family, K, eps, delta, L = NULL)
}
\arguments{
\item{x}{A matrix or data frame of predictor variables. Column names will be used
as variable names; if missing, variables will be named X1, X2, etc.}

\item{y}{A vector of response values. Should be numeric for linear regression
or binary (0/1) for logistic regression.}

\item{family}{A character string specifying the model family. Must be either
"gaussian" for linear regression or "binomial" for logistic regression.}

\item{K}{An integer specifying the maximum number of steps (variables to add).
The algorithm will terminate early if no models meet the improvement criteria.}

\item{eps}{A numeric threshold for minimum AIC improvement required to add a
variable. Variables must improve parent model AIC by at least this amount.}

\item{delta}{A numeric threshold for competitiveness. Child models are kept if
their AIC is within delta of the best child model at that step.}

\item{L}{An integer specifying the maximum number of models to retain at each
step (frontier size). If NULL, all qualifying models are kept. Use this to
control computational cost for high-dimensional problems.}
}
\value{
A list with class "path_forest" containing:
\describe{
\item{frontiers}{A list of data frames, one per step, containing the models
at each frontier. Each data frame has columns: step, variables (list column),
aic, model_id, and parent_id.}
\item{aic_by_model}{A named list storing AIC values for all explored models,
indexed by model_id to avoid redundant computation.}
\item{meta}{A list of metadata including family, K, eps, delta, L, number of
steps completed, number of variables, variable names, and total models explored.}
}
}
\description{
Constructs a forest of model selection paths by iteratively adding variables
to models based on AIC improvement criteria. At each step, the algorithm explores
multiple promising paths simultaneously, keeping models that improve sufficiently
and are competitive with the best model at that step.
}
\details{
The algorithm starts with an empty model (intercept only) and iteratively builds
a forest of model selection paths. At each step k:
\enumerate{
\item For each model in the current frontier, evaluate adding each available variable
\item Keep child models that: (a) improve parent AIC by at least eps, AND
(b) are within delta AIC of the best child at this step
\item Stop if no models meet criteria or K steps are reached
}

The multi-path approach explores multiple promising directions simultaneously,
unlike traditional stepwise selection which follows a single path.
}
\examples{
# Linear regression example
set.seed(123)
n <- 100
x <- matrix(rnorm(n * 5), ncol = 5)
colnames(x) <- paste0("Var", 1:5)
y <- x[,1] + 2*x[,2] + rnorm(n)

paths <- build_paths(x, y, family = "gaussian", K = 5, eps = 0.5, delta = 2)

# Logistic regression example
y_binary <- rbinom(n, 1, plogis(x[,1] + x[,2]))
paths_logit <- build_paths(x, y_binary, family = "binomial", 
                          K = 5, eps = 0.5, delta = 2, L = 10)
}
\seealso{
\code{\link{stability}} for analyzing path stability,
\code{\link{plausible_models}} for choosing final models
}

