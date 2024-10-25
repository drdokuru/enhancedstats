\name{enhancedSummary}
\alias{enhancedSummary}
\title{Enhanced Summary Function for Linear and Mixed-Effects Models}
\description{
  This function provides a comprehensive summary of linear and mixed-effects models,
  including model performance metrics, ANOVA table, coefficient estimates, and multicollinearity diagnostics.
}
\usage{
enhancedSummary(model, correlation = FALSE, symbolic.cor = FALSE, ...)
}
\arguments{
  \item{model}{A fitted model of class 'lm' (linear model) or 'merMod' (mixed-effects model).}
  \item{correlation}{Logical; if TRUE, calculates and displays correlation of coefficients (default = FALSE).}
  \item{symbolic.cor}{Logical; if TRUE, displays a symbolic representation of the correlation matrix (default = FALSE).}
  \item{...}{Additional arguments passed to specific model functions.}
}
\value{
  Prints a detailed summary of the model including ANOVA table, model performance,
  coefficient estimates, and correlation of coefficients (if requested).
}
\examples{
data(mtcars)
lm_model <- lm(mpg ~ wt + hp, data = mtcars)
enhancedSummary(lm_model)

# For mixed-effects models
# library(lme4)
# mixed_model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
# enhancedSummary(mixed_model)
}
