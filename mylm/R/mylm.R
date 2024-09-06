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


  TSS <- sum((y-mean(y))^2)

  # and store the results in the list est
  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula
  est$coeff <- coeff_list
  est$rank <- length(colnames(X))
  est$fitted_values <- fitted_values
  est$residuals <- residuals

  est$dof_residuals <- nrow(X) - length(colnames(X))
  est$data_matrix <- X
  est$TSS <- TSS


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

  #
  X <- object$data_matrix
  RSS <- sum(object$residuals^2)
  sigma2 <- RSS/object$dof_residuals
  XtX_inv <- solve(t(X)%*%X)
  cov_matrix <- sigma2*XtX_inv
  stderr <- sqrt(diag(cov_matrix))


  # z and p values
  z <- as.numeric(object$coeff) / as.numeric(stderr)
  p <- 2 * (1 - pnorm(abs(z)))

  # significance levels
  sig_level <- list()
  for (value in p) {
    # Determine the significance level and append to the list
    if (value < 0.001) {
      sig_level[[length(sig_level) + 1]] <- '***'
    } else if (value < 0.01) {
      sig_level[[length(sig_level) + 1]] <- '**'
    } else if (value < 0.05) {
      sig_level[[length(sig_level) + 1]] <- '*'
    } else if (value < 0.1) {
      sig_level[[length(sig_level) + 1]] <- '.'
    } else {
      sig_level[[length(sig_level) + 1]] <- ' '
    }
  }


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
  formatted_residuals <- format(summary_residuals, digits = 4, nsmall = 3,justify = "right", trim = TRUE)
  max_width <- max(nchar(formatted_residuals))
  formatted_residuals <- format(summary_residuals, digits = 4, nsmall = 3, justify = "right", trim = TRUE)

  # printing
  cat('\nResiduals:\n')

  cat(str_pad("Min", max_width+2, side = 'right'),
      str_pad("1Q", max_width+2, side = 'right'),
      str_pad("Median", max_width+2, side = 'right'),
      str_pad("3Q", max_width+2, side = 'right'),
      str_pad("Max", max_width+2, side = 'right'), "\n")
  cat(str_pad(formatted_residuals[1], max_width+2, side = 'right'),
      str_pad(formatted_residuals[2], max_width+2, side = 'right'),
      str_pad(formatted_residuals[3], max_width+2, side = 'right'),
      str_pad(formatted_residuals[4], max_width+2, side = 'right'),
      str_pad(formatted_residuals[5], max_width+2, side = 'right'), "\n")




  cat('\nCoefficients:\n')
  max_name = max(nchar(names(object$coeff)))
  formatted_coeff<- format(object$coeff, digits = 4, nsmall = 6,justify = "right", trim = TRUE)
  max_width <- max(nchar(formatted_coeff))
  formatted_coeff <- format(object$coeff, digits = 4, nsmall = 6, justify = "right", trim = TRUE)


  cat(strrep(" ", max_name+2),
      str_pad('Estimate', max_width+3, 'right'),
      str_pad('Std. Error', max_width+3, 'right'),
      str_pad("z value", max_width+3, 'right'),
      str_pad( "Pr(>|z|)", max_width+3, 'right'), '\n')
  i <- 1
  for (name in names(object$coeff)) {
    cat(str_pad(name, max_name+3, 'right'))
    cat(
        str_pad(formatted_coeff[[name]], max_width+3, 'right'),
        str_pad(format(stderr[i], digits = 4, nsmall = 6, justify = "right", trim = TRUE), max_width+3, 'right'),
        str_pad(format(z[i], digits = 4, nsmall = 6, justify = "right", trim = TRUE), max_width+3, 'right'),
        str_pad(paste(format(p[i], digits = 4, nsmall = 6, justify = "right", trim = TRUE), sig_level[i]), max_width+3, 'right'),
        '\n')
    i <- i+1
  }

  R_sqrd <- 1-(RSS/object$TSS)
  F_stat <- ((1/object$rank-1)*(object$TSS-RSS))/(RSS/object$dof_residuals )
  #p_F <-

  cat('\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n')

  cat('Residual standard error: ', format(sqrt(sigma2), digits=4), 'on ', object$dof_residuals,'degrees of freedom', '\n') # values missing
  cat('Multiple R-squared: ', format(R_sqrd, digits=4), ', ', 'Adjusted R-squared: ', format(1-((1-R_sqrd)*(nrow(X)-1)/object$dof_residuals), digits=4), "\n")
  cat('F-statistic: ', F_stat, 'on',object$rank-1, 'and',object$dof_residuals ,'DF, p-value:' )

}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"

  library(ggplot2)
  # ggplot requires that the data is in a data.frame, this must be done here
  data_plot <- data.frame(Fitted = object$fitted_values, Residuals=object$residuals)

    ggplot(data_plot, aes(x=Fitted, y=Residuals)) + geom_point(shape = 1, size = 2) +
      ggtitle('Residuals vs Fitted') +
      labs(x = 'Fitted Values', y = 'Residuals') +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      theme_minimal()

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
