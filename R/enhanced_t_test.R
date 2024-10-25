#' Enhanced T-Test
#'
#' This function performs enhanced t-tests (one-sample, paired, two-sample) and reports effect sizes.
#'
#' @param x Numeric vector of data values.
#' @param y Numeric vector of data values (optional for one-sample t-test).
#' @param alpha Significance level (default = 0.05).
#' @param paired Logical; if TRUE, performs a paired t-test (default = FALSE).
#' @param mu Population mean (default = NULL).
#' @param var.equal Logical; if TRUE, assumes equal variances for two-sample t-tests (default = TRUE).
#'
#' @return Prints the results of the t-test, including statistics and effect sizes.
#' @examples
#' enhanced_t_test(c(5, 6, 7), mu = 6)
#' enhanced_t_test(c(5, 6, 7), c(6, 7, 8), paired = TRUE)
#' enhanced_t_test(c(5, 6, 7), c(8, 9, 10), var.equal = FALSE)
enhanced_t_test <- function(x, y = NULL, alpha = 0.05, paired = FALSE, mu = NULL, var.equal = TRUE) {

  # Interpret Cohen's d
  interpret_cohen_d <- function(d) {
    if (is.na(d)) return("Effect size not applicable")
    if (abs(d) < 0.2) {
      return("Small effect")
    } else if (abs(d) < 0.5) {
      return("Medium effect")
    } else if (abs(d) < 0.8) {
      return("Large effect")
    } else {
      return("Very large effect")
    }
  }

  # Determine the type of t-test to perform
  if (is.null(y)) {
    # One-Sample T-Test
    t_test <- t.test(x, mu = mu, conf.level = 1 - alpha)
    null_hypothesis <- sprintf("H0: The mean of the sample is equal to %.2f", mu)

    # Statistics for one-sample
    mean_x <- mean(x)
    sd_x <- sd(x)
    se_x <- sd_x / sqrt(length(x))

    # Raw Effect Size and Cohen's d
    raw_effect_size <- mean_x - mu
    cohen_d <- raw_effect_size / sd_x

    # Cohen's d Interpretation
    cohen_d_interpretation <- interpret_cohen_d(cohen_d)

    # Output
    cat("\033[1mOne-Sample T-Test Results:\033[0m\n")
    cat("-----------------------------------------\n")
    cat("\033[1mNull Hypothesis:\033[0m\n")
    cat(null_hypothesis, "\n\n")
    cat("\033[1mGroup Statistics:\033[0m\n")
    cat(sprintf("Mean = %.3f, SD = %.3f, SE = %.3f\n", mean_x, sd_x, se_x))
    cat("-----------------------------------------\n")
    cat("\033[1mEffect Size:\033[0m\n")
    cat(sprintf("Raw Effect Size (Mean Difference) = %.3f\n", raw_effect_size))
    cat(sprintf("Cohen's d = %.3f (%s)\n", cohen_d, cohen_d_interpretation))
    cat("-----------------------------------------\n")
    cat("\033[1m95%% Confidence Interval:\033[0m\n")
    cat(sprintf("[%.3f, %.3f]\n", t_test$conf.int[1], t_test$conf.int[2]))
    cat("-----------------------------------------\n")
    cat("\033[1mTest Statistics:\033[0m\n")
    cat(sprintf("t(%.1f) = %.3f, p-value = %.4f\n", t_test$parameter, t_test$statistic, t_test$p.value))

  } else if (paired) {
    # Paired T-Test
    t_test <- t.test(x, y, paired = TRUE, conf.level = 1 - alpha)
    null_hypothesis <- "H0: The means of the paired groups are equal (Mx = My)"

    # Statistics for paired
    mean_x <- mean(x)
    mean_y <- mean(y)
    sd_x <- sd(x)
    sd_y <- sd(y)

    # Raw Effect Size and Cohen's d
    raw_effect_size <- mean_x - mean_y
    cohen_d <- raw_effect_size / sd(x)  # Using sd(x) since pairs are related

    # Cohen's d Interpretation
    cohen_d_interpretation <- interpret_cohen_d(cohen_d)

    # Output
    cat("\033[1mPaired T-Test Results:\033[0m\n")
    cat("-----------------------------------------\n")
    cat("\033[1mNull Hypothesis:\033[0m\n")
    cat(null_hypothesis, "\n\n")
    cat("\033[1mGroup Statistics:\033[0m\n")
    cat(sprintf("Group 1: Mean = %.3f, SD = %.3f\n", mean_x, sd_x))
    cat(sprintf("Group 2: Mean = %.3f, SD = %.3f\n", mean_y, sd_y))
    cat("-----------------------------------------\n")
    cat("\033[1mEffect Size:\033[0m\n")
    cat(sprintf("Raw Effect Size (Mean Difference) = %.3f\n", raw_effect_size))
    cat(sprintf("Cohen's d = %.3f (%s)\n", cohen_d, cohen_d_interpretation))
    cat("-----------------------------------------\n")
    cat("\033[1m95%% Confidence Interval:\033[0m\n")
    cat(sprintf("[%.3f, %.3f]\n", t_test$conf.int[1], t_test$conf.int[2]))
    cat("-----------------------------------------\n")
    cat("\033[1mTest Statistics:\033[0m\n")
    cat(sprintf("t(%.1f) = %.3f, p-value = %.4f\n", t_test$parameter, t_test$statistic, t_test$p.value))

  } else {
    # Two-Sample T-Test
    t_test <- t.test(x, y, var.equal = var.equal, conf.level = 1 - alpha)
    null_hypothesis <- "H0: The means of the two groups are equal (Mx = My)"

    # Statistics for two-sample
    mean_x <- mean(x)
    sd_x <- sd(x)
    se_x <- sd_x / sqrt(length(x))

    mean_y <- mean(y)
    sd_y <- sd(y)
    se_y <- sd_y / sqrt(length(y))

    # Pooled Standard Deviation and SE of Difference
    n_x <- length(x)
    n_y <- length(y)

    if (var.equal) {
      # Pooled Standard Deviation when variances are assumed equal
      pooled_sd <- sqrt(((n_x - 1) * sd_x^2 + (n_y - 1) * sd_y^2) / (n_x + n_y - 2))
      se_diff <- sqrt((sd_x^2 / n_x) + (sd_y^2 / n_y))
    } else {
      # Standard Deviations for separate variances
      se_diff <- sqrt((sd_x^2 / n_x) + (sd_y^2 / n_y))
      pooled_sd <- NA  # Not applicable when variances are unequal
    }

    # Raw Effect Size and Cohen's d
    raw_effect_size <- mean_x - mean_y
    cohen_d <- if (!is.na(pooled_sd)) raw_effect_size / pooled_sd else NA

    # Cohen's d Interpretation
    cohen_d_interpretation <- interpret_cohen_d(cohen_d)

    # Output
    cat("\033[1mTwo-Sample T-Test Results:\033[0m\n")
    cat("-----------------------------------------\n")
    cat("\033[1mNull Hypothesis:\033[0m\n")
    cat(null_hypothesis, "\n\n")
    cat("\033[1mGroup Statistics:\033[0m\n")
    cat(sprintf("Group 1: Mean = %.3f, SD = %.3f, SE = %.3f\n", mean_x, sd_x, se_x))
    cat(sprintf("Group 2: Mean = %.3f, SD = %.3f, SE = %.3f\n", mean_y, sd_y, se_y))
    cat("-----------------------------------------\n")
    cat("\033[1mEffect Size:\033[0m\n")
    cat(sprintf("Raw Effect Size (Mean Difference) = %.3f\n", raw_effect_size))
    cat(sprintf("Cohen's d = %.3f (%s)\n", cohen_d, cohen_d_interpretation))
    cat("-----------------------------------------\n")

    if (!is.null(pooled_sd)) {
      cat("\033[1mPooled Standard Deviation:\033[0m\n")
      cat(sprintf("Pooled SD = %.3f\n", pooled_sd))
      cat("-----------------------------------------\n")
    }

    cat("\033[1mStandard Error of the Difference:\033[0m\n")
    cat(sprintf("SE of Difference in Means = %.3f\n", se_diff))
    cat("-----------------------------------------\n")
    cat("\033[1m95%% Confidence Interval:\033[0m\n")
    cat(sprintf("[%.3f, %.3f]\n", t_test$conf.int[1], t_test$conf.int[2]))
    cat("-----------------------------------------\n")
    cat("\033[1mTest Statistics:\033[0m\n")
    cat(sprintf("t(%.1f) = %.3f, p-value = %.4f\n", t_test$parameter, t_test$statistic, t_test$p.value))
  }
}

# Example usage (uncomment to run)
# enhanced_t_test(c(5, 6, 7), mu = 6)
# enhanced_t_test(c(5, 6, 7), c(6, 7, 8), paired = TRUE)
# enhanced_t_test(c(5, 6, 7), c(8, 9, 10), var.equal = FALSE)
