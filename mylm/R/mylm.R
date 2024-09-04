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
  coeff_list <- as.list(coeff[,1])

  # Assign names to coefficients
  names(coeff_list) <- colnames(X)

  # fitted values
  fitted_values <- X %*% coeff
  #residuals
  residuals <- y - fitted_values


  # and store the results in the list est
  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula
  est$coeff <- coeff_list
  est$rank <- length(colnames(X))
  est$fitted_values <- fitted_values
  est$residuals <- residuals
  est$df_residuals <- nrow(X) - length(colnames(X))



  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  variable_names = all.vars(object$formula)

  cat('Call:\n')
  print(object$call)
  cat('\nCoefficients:\n')
  for (name in names(object$coeff)) {
    cat(name, ': ')
    cat(format(object$coeff[[name]], digits = 4), '\n')
  }
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format

  ## Call
  cat('Call:\n')
  print(object$call)

  ## Residuals
  # set up values
  summary_residuals <- c(
    Min = min(object$residuals),
    Q1 = quantile(object$residuals, 0.25),
    Median = median(object$residuals),
    Q3 = quantile(object$residuals, 0.75),
    Max = max(object$residuals)
  )
  formatted_residuals <- format(summary_residuals, digits = 4, nsmall = 3, justify = "right")
  max_width <- max(nchar(formatted_residuals))
  formatted_residuals <- format(summary_residuals, digits = 4, nsmall = 3, justify = "right")

  # printing
  cat('\nResiduals:\n')
  cat("Min", strrep(" ", max_width - nchar("Min")+2), "1Q",
      strrep(" ", max_width - nchar("1Q")+2), "Median",
      strrep(" ", max_width - nchar("Median")+2), "3Q",
      strrep(" ", max_width - nchar("3Q")+2), "Max\n")
  cat(formatted_residuals[1], " ", formatted_residuals[2], " ",
      formatted_residuals[3], " ", formatted_residuals[4], " ",
      formatted_residuals[5], "\n")



  cat('\nCoefficients:\n')


  print('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1')

  cat('Residual standard error: ', 'on ', 'degrees of freedom', '\n') # values missing
  cat('Multiple R-squared: ', ', ', 'Adjusted R-squared: ' )
  cat(' F-statistic: ', 'on', 'and', 'DF, p-value:', )

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
