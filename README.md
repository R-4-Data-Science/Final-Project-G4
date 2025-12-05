Example Linear Regression Model

```{r}

#install package and load library
remotes::install_github("R-4-Data-Science/Final-Project-G4")
library("MultiPathSelection")


##Linear regression
set.seed(1)
n <- 120; p <- 8
X <- matrix(rnorm(n*p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p-5))
y <- X %*% beta + rnorm(n, sd = 1)
colnames(X) <- paste0("x", 1:p)
df <- as.data.frame(cbind(y, X))

# Fit a few single-step candidates to illustrate AIC comparisons:
aic_empty <- AIC(lm(y ~ 1, data = df))
cand <- lapply(1:p, function(j) lm(as.formula(paste0("y ~ ", colnames(X)[j])), data = df))
aics <- sapply(cand, AIC)
data.frame(variable = colnames(X), AIC = aics)[order(aics), ][1:5, ]

#use build_paths function to create a forest of possible models
forest <- build_paths(
  x = X,
  y = as.numeric(y),
  family = "gaussian",
  K = min(ncol(X), 10),
  eps = 1e-6,
  delta = 1,
  L = 50
)
# use stability() function to find reappearing predictors to determine stability for each model
stab <- stability(
  x = X,
  y = as.numeric(y),
  B = 50,
  resample = "bootstrap",
  family = "gaussian",
  K = 10,
  eps = 1e-6,
  delta = 1,
  L = 50
)

#use plausible_models function to sort out the less stable models and return a matrix with stats of the plausible models
# used higher delta and lower tau values becuase the program was too strict and was returning no functions
plaus <- plausible_models(
  forest,
  pi = stab$pi,
  Delta = 6,
  tau = 0.3
)
plaus
# x1+x2+x5+x7 and x1+x2+x5 are the plausible models for our arbitrary linear regression.
```
Example Logistic Regression
```{r}
#install package and load library
remotes::install_github("R-4-Data-Science/Final-Project-G4")
library("MultiPathSelection")

#perform logistic regression on an arbitrary dataset
set.seed(2)
n <- 200; p <- 6
Xb <- matrix(rnorm(n*p), n, p)
linpred <- 1.2*Xb[,1] - 1*Xb[,2] + 0.8*Xb[,5]
prob <- 1 / (1 + exp(-linpred))
ybin <- rbinom(n, 1, prob)
colnames(Xb) <- paste0("x", 1:p)
dfb <- as.data.frame(cbind(y = ybin, Xb))

#
fit0 <- glm(y ~ 1, family = binomial(), data = dfb)
aic0 <- AIC(fit0)
fits1 <- lapply(1:p, function(j) glm(as.formula(paste0("y ~ ", colnames(Xb)[j])),
                                     family = binomial(), data = dfb))
aics1 <- sapply(fits1, AIC)
head(data.frame(variable = colnames(Xb), AIC = aics1)[order(aics1), ], 5)

# use build_paths() function to create a forest varible with all possible models to use
forest <- build_paths(
  x = Xb,
  y = ybin,
  family = "binomial",
  K = min(ncol(Xb), 5),
  eps = 1e-6,
  delta = 1,
  L = 50
)

#use stability() function to record stability metrics for each model
stab <- stability(
  x = Xb,
  y = ybin,
  B = 50,
  resample = "bootstrap",
  family = "binomial",
  K = 5,
  eps = 1e-6,
  delta = 1,
  L = 50
)

#filter out the low-stability models and return the best fitted ones
plaus <- plausible_models(
  forest,
  pi = stab$pi,
  Delta = 2,
  tau = 0.3
)

print(plaus)
# for the logistic data, we find that x1+x2+x3+x5 and x1+x2+x5 are plausible models for the logistic data
```
We were left with 2 models each for the linear and logistic regressions.
