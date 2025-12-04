# takes list of variable input and returns a list with model and AIC value
fit_model <- function(variables, x, y, family) {
  #model with one or more predictors
  formula_str <- paste("y ~", paste(variables, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  if (family == "gaussian") {
    model <- lm(formula_obj, data = x)
  } else {
    model <- glm(formula_obj, data = x, family = binomial)
  }

  return(list(model = model, aic = AIC(model)))
}

###
# builds empty model and returns a list with model and AIC value
fit_empty_model <- function(y, family) {
  if (family == "gaussian") {
    model <- lm(y ~ 1)
  } else {
    model <- glm(y ~ 1, family = binomial)
  }

  return(list(model = model, aic = AIC(model)))
}

###
# create unique model identifier
model_id <- function(variables) {
  if (length(variables) == 0) return("intercept_only")
  paste(sort(variables), collapse = "+")
}
