# Select Build, Build and reload to build and lode into the R-session.

library(stringr)

mylm <- function(formula, data = list(), contrasts = NULL, ...){
  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")


  # Add code here to calculate coefficients, residuals, fitted values, etc...
  # coefficients
  coeff <- solve(t(X) %*% X) %*% t(X) %*% y
  # fitted values
  fitted_values <- X %*% coeff
  #residuals
  residuals <- y - fitted_values


  # and store the results in the list est
  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula
  est$coeff <- coeff

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  variable_names = all.vars(object$formula)
  #max_length <- max(lapply(variable_names, length))
  print(variable_names)

  cat('Call:\n')
  print(object$call)
  cat('Coefficients:\n')
  cat('Intercept : ')
  cat(format(object$coeff[1], digits = 4), '\n')
  for(i in 2:length(variable_names)){
    cat(str_pad(variable_names[i], width = max_length + 2, side = "right", pad = " "),': ')
    cat(format(object$coeff[i], digits = 4), '\n')

  }
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Summary of object\n')
}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"

  library(ggplot2)
  # ggplot requires that the data is in a data.frame, this must be done here
  ggplot() + geom_point()

  # if you want the plot to look nice, you can e.g. use "labs" to add labels, and add colors in the geom_point-function

}

anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- lm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
  }

  return(model)

}
