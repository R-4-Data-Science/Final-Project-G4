---
output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  message = FALSE,
  warning = FALSE
)
```

# MultiPathSelection

<!-- badges: start -->
<!-- badges: end -->

MultiPathSelection implements a multi-path selection algorithm for building robust linear and logistic regression models. Unlike traditional stepwise selection that follows a single path, this package explores multiple promising model-building paths simultaneously, using stability analysis to identify reliable variable selections.

## Key Features

- **Multi-path exploration**: Builds multiple model paths simultaneously instead of a single greedy path
- **Stability-based selection**: Uses resampling to identify variables that are consistently selected
- **AIC-guided search**: Balances model fit and complexity using information criteria
- **Works for both**: Linear regression (Gaussian) and logistic regression (binomial)

## Installation

Install the development version from GitHub:
```{r, eval = FALSE}
# install.packages("remotes")
remotes::install_github("R-4-Data-Science/Final-Project-G4")
```

## Quick Start

The typical workflow involves three main steps:

1. **Build paths** - Explore multiple model-building trajectories
2. **Compute stability** - Identify consistently selected variables via resampling
3. **Select plausible models** - Filter to the most stable, competitive models

## Example: Linear Regression
```{r linear-example}
library(MultiPathSelection)

# Simulate data with true model: y = 2*x1 - 1.5*x2 + x5 + noise
set.seed(1)
n <- 120
p <- 8
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
y <- X %*% beta + rnorm(n, sd = 1)

# Step 1: Build path forest
forest <- build_paths(
  x = X,
  y = as.numeric(y),
  family = "gaussian",
  K = 8,           # Maximum steps
  eps = 1e-6,      # Minimum AIC improvement
  delta = 1,       # Competitiveness threshold
  L = 50           # Max models per step
)

# Step 2: Compute stability via bootstrap
stab <- stability(
  x = X,
  y = as.numeric(y),
  B = 50,              # 50 bootstrap samples
  resample = "bootstrap",
  family = "gaussian",
  K = 8,
  eps = 1e-6,
  delta = 1,
  L = 50
)

# View stability scores
print(sort(stab$pi, decreasing = TRUE))

# Step 3: Select plausible models
plaus <- plausible_models(
  forest,
  pi = stab$pi,
  Delta = 4,       # AIC competitiveness
  tau = 0.5        # Minimum avg stability
)

print(plaus)
```

The algorithm correctly identifies the important variables (x1, x2, x5) from the true model!

## Example: Logistic Regression
```{r logistic-example}
# Simulate binary response data
set.seed(2)
n <- 200
p <- 6
X_log <- matrix(rnorm(n * p), n, p)
colnames(X_log) <- paste0("x", 1:p)

# True model: logit(P(y=1)) = 1.2*x1 - x2 + 0.8*x5
linpred <- 1.2 * X_log[, 1] - 1 * X_log[, 2] + 0.8 * X_log[, 5]
prob <- 1 / (1 + exp(-linpred))
y_bin <- rbinom(n, 1, prob)

# Build paths for logistic regression
forest_log <- build_paths(
  x = X_log,
  y = y_bin,
  family = "binomial",  # Logistic regression
  K = 6,
  eps = 1e-6,
  delta = 1,
  L = 50
)

# Compute stability
stab_log <- stability(
  x = X_log,
  y = y_bin,
  B = 50,
  resample = "bootstrap",
  family = "binomial",
  K = 6,
  eps = 1e-6,
  delta = 1,
  L = 50
)

print(sort(stab_log$pi, decreasing = TRUE))

# Select plausible models
plaus_log <- plausible_models(
  forest_log,
  pi = stab_log$pi,
  Delta = 2,
  tau = 0.5
)

print(plaus_log)
```

## Confusion Metrics for Classification

For logistic regression, evaluate model performance:
```{r confusion}
# Fit the best model
best_vars <- plaus_log$variables[[1]]
formula <- as.formula(paste("y_bin ~", paste(best_vars, collapse = " + ")))
final_model <- glm(formula, data = as.data.frame(X_log), family = binomial())

# Get predictions
y_prob <- predict(final_model, type = "response")

# Compute confusion matrix and metrics
metrics <- confusion_metrics(y_bin, y_prob, cutoff = 0.5)
print(metrics$confusion_matrix)
print(metrics$metrics)
```

## Main Functions

- `build_paths()` - Construct forest of model selection paths
- `stability()` - Compute variable stability scores via resampling  
- `plausible_models()` - Filter to stable, competitive models
- `confusion_metrics()` - Evaluate binary classification performance

## Getting Help

View detailed documentation and examples:
```{r, eval = FALSE}
?build_paths
?stability
?plausible_models
browseVignettes("MultiPathSelection")
```

## License

MIT License

## Authors

Jonah Kennedy (jak0110@auburn.edu)
