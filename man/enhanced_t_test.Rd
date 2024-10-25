\name{enhanced_t_test}
\alias{enhanced_t_test}
\title{Enhanced T-Test Function}
\description{
  This function performs enhanced t-tests (one-sample, paired, two-sample) and reports effect sizes.
}
\usage{
enhanced_t_test(x, y = NULL, alpha = 0.05, paired = FALSE, mu = NULL, var.equal = TRUE)
}
\arguments{
  \item{x}{A numeric vector of data values.}
  \item{y}{A numeric vector of data values (optional for one-sample t-test).}
  \item{alpha}{A numeric value specifying the significance level (default is 0.05).}
  \item{paired}{Logical; if TRUE, performs a paired t-test (default is FALSE).}
  \item{mu}{A numeric value representing the population mean (default is NULL).}
  \item{var.equal}{Logical; if TRUE, assumes equal variances for two-sample t-tests (default is TRUE).}
}
\value{
  Prints the results of the t-test, including test statistics and effect sizes.
}
\examples{
# One-sample t-test
enhanced_t_test(c(5, 6, 7), mu = 6)

# Two-sample t-test
enhanced_t_test(c(5, 6, 7), c(7, 8, 9), var.equal = TRUE)

# Paired t-test
enhanced_t_test(c(5, 6, 7), c(6, 7, 8), paired = TRUE)
}
