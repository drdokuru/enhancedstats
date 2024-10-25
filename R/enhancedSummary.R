#' Enhanced Summary Function for Linear and Mixed-Effects Models
#'
#' This function provides a comprehensive summary of linear and mixed-effects models,
#' including model performance metrics, ANOVA table, coefficient estimates, and multicollinearity diagnostics.
#'
#' @param model A fitted model of class 'lm' (linear model) or 'merMod' (mixed-effects model).
#' @param correlation Logical; if TRUE, calculates and displays correlation of coefficients (default = FALSE).
#' @param symbolic.cor Logical; if TRUE, displays a symbolic representation of the correlation matrix (default = FALSE).
#' @param ... Additional arguments passed to specific model functions.
#'
#' @return Prints a detailed summary of the model including ANOVA table, model performance,
#' coefficient estimates, and correlation of coefficients (if requested).
#' @examples
#' data(mtcars)
#' lm_model <- lm(mpg ~ wt + hp, data = mtcars)
#' enhancedSummary(lm_model)
#'
#' # For mixed-effects models
#' # library(lme4)
#' # mixed_model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' # enhancedSummary(mixed_model)
enhancedSummary <- function(model, correlation = FALSE, symbolic.cor = FALSE, ...) {

  # Check if the 'car' package is available, install it if not
  if (!requireNamespace("car", quietly = TRUE)) {
    install.packages('car')
  }
  library(car)

  # Check if the model is of class 'lm' (linear model)
  if (inherits(model, "lm")) {
    summary_model <- summary(model)  # Retrieve summary of the model
    m <- model  # Assign the model to 'm'
    p <- m$rank  # Extract rank of the model matrix (design matrix)
    rdf <- m$df.residual  # Extract residual degrees of freedom
    Qr <- m$qr  # QR decomposition of the model matrix
    n <- NROW(Qr$qr)  # Total number of observations

    # Check if there are predictors in the model
    if (p == 0) {
      stop("No predictors in the model")
    }

    # Calculations
    r <- m$residuals  # Residuals
    f <- m$fitted.values  # Fitted values
    w <- m$weights  # Observation weights

    # If no weights are specified, use equal weights
    if (is.null(w)) {
      mss <- if (attr(m$terms, "intercept")) sum((f - mean(f))^2) else sum(f^2)  # Model sum of squares
      rss <- sum(r^2)  # Residual sum of squares
    } else {
      mss <- if (attr(m$terms, "intercept")) {
        m <- sum(w * f/sum(w))
        sum(w * (f - m)^2)
      } else sum(w * f^2)
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }

    resvar <- rss / rdf  # Residual variance
    r.squared <- mss / (mss + rss)  # Compute R-squared
    adj.r.squared <- 1 - (1 - r.squared) * ((n - p) / rdf)  # Compute adjusted R-squared
    rmse <- sqrt(resvar)  # Calculate root mean squared error (RMSE)

    # ANOVA table
    df_model <- p - 1
    df_residual <- rdf
    df_total <- n - 1
    SS_model <- mss
    SS_residual <- rss
    SS_total <- mss + rss
    MS_model <- SS_model / df_model
    MS_residual <- resvar
    F_value <- MS_model / MS_residual
    p_value <- pf(F_value, df_model, df_residual, lower.tail = FALSE)

    # Create ANOVA table
    anova_table <- data.frame(
      Source = c("Model", "Residuals", "Total"),
      SS = c(SS_model, SS_residual, SS_total),
      df = c(df_model, df_residual, df_total),
      MS = c(MS_model, MS_residual, NA),
      F = c(F_value, NA, NA),
      p = c(p_value, NA, NA)
    )

    # Coefficients table
    coef_table <- summary_model$coefficients
    coef_table <- cbind(coef_table, confint(model))
    colnames(coef_table)[5:6] <- c("CI(2.5)", "CI(97.5)")
    colnames(coef_table)[4] <- "p"  # Renaming the p-value column to "p"

    # VIF for multicollinearity
    vif_values <- vif(model)
    tol <- 1 / vif_values
    coef_table <- cbind(coef_table, tol)
    colnames(coef_table)[7] <- "Tol"

    # Output
    print(m$call)
    cat("\nANOVA Table\n")
    print(anova_table, row.names = FALSE)
    cat("\nModel Performance\n")
    cat(sprintf("RMSE: %.3f\n", rmse))
    cat(sprintf("Adjusted R-squared: %.3f\n", adj.r.squared))
    cat("\nCoefficients\n")
    print(round(coef_table, 3))

    if (correlation) {
      cor_matrix <- cor(m$model)
      cat("\nCorrelation of Coefficients:\n")
      print(round(cor_matrix, 3))
      if (symbolic.cor) {
        print(symnum(cor_matrix))
      }
    }

  } else if (inherits(model, "merMod")) {  # Check if the model is of class 'merMod' (mixed-effects model)

    if (!requireNamespace("lme4", quietly = TRUE)) {
      install.packages('lme4')
    }
    library(lme4)
    if (!requireNamespace("MuMIn", quietly = TRUE)) {
      install.packages('MuMIn')
    }
    if (!requireNamespace("lmerTest", quietly = TRUE)) {
      install.packages('lmerTest')
    }
    library(lmerTest)

    summary_model <- summary(as(model, "merModLmerTest"))  # Retrieve summary of the model
    fixed_effects <- summary_model$coefficients  # Extract fixed effects coefficients
    random_effects <- summary_model$varcor  # Extract random effects variance components
    n <- length(model@resp$y)  # Total number of observations
    p <- ncol(fixed_effects)  # Number of fixed effects parameters
    rdf <- n - p  # Residual degrees of freedom

    # Calculate RMSE
    residuals <- resid(model)
    rmse <- sqrt(mean(residuals^2))

    # Calculate R-squared and adjusted R-squared
    r.squared <- r.squaredGLMM(model)
    marginal_r2 <- r.squared[1, "R2m"]
    conditional_r2 <- r.squared[1, "R2c"]

    # Coefficients table with confidence intervals and p-values
    confint_vals <- as.data.frame(confint(model))
    ci_fixed <- confint_vals[rownames(fixed_effects), ]

    # Check if fixed_effects has at least four columns
    if (ncol(fixed_effects) >= 4) {
      coef_table <- cbind(
        Estimate = fixed_effects[, 1],
        Std.Error = fixed_effects[, 2],
        t.value = fixed_effects[, 4],
        p = fixed_effects[, 5],  # Access the p-values
        'CI(2.5)' = ci_fixed[, 1],
        'CI(97.5)' = ci_fixed[, 2]
      )
    } else {
      coef_table <- cbind(
        Estimate = fixed_effects[, 1],
        Std.Error = fixed_effects[, 2],
        t.value = fixed_effects[, 3],
        p = NA,  # Set p-values to NA if unavailable
        'CI(2.5)' = ci_fixed[, 1],
        'CI(97.5)' = ci_fixed[, 2]
      )
    }

    # Output
    print(model@call)
    cat("\nRandom Effects:\n")
    print(random_effects)
    cat("\nModel Performance\n")
    cat(sprintf("RMSE: %.3f\n", rmse))
    cat(sprintf("Marginal R-squared: %.3f\n", marginal_r2))
    cat(sprintf("Conditional R-squared: %.3f\n", conditional_r2))
    cat("\nFixed Effects Coefficients\n")
    print(coef_table)

    if (correlation) {
      cor_matrix <- cov2cor(vcov(model))
      cat("\nCorrelation of Fixed Effects Coefficients:\n")
      print(round(cor_matrix, 3))
      if (symbolic.cor) {
        print(symnum(cor_matrix))
      }
    }

  } else {
    stop("The model must be of class 'lm' or 'merMod'")
  }
}

# Example usage (uncomment to run)
# library(lme4)
# mixed_model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
# enhancedSummary(mixed_model)
