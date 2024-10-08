% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preGroup.r
\name{coef.preGroup}
\alias{coef.preGroup}
\title{Extract Grouped Coefficients from preGroup Model}
\usage{
\method{coef}{preGroup}(object, ...)
}
\arguments{
\item{object}{A model object, which must be of class \code{preGroup}.}

\item{...}{Additional arguments passed to or from other methods.}
}
\value{
A data frame of combined coefficients and descriptions, sorted with respect to origin.
}
\description{
This function calculates the aggregated coefficients for a \code{preGroup} model,
combining intercepts and coefficients from different hierarchical levels.
}
\details{
The function first verifies the class of the model, then extracts and processes
coefficients from multiple levels of the model to compute the grand intercept and other
subgroup coefficients. It finally merges these results with the rule descriptions from
the \code{pre_fit} component of the model, handles missing descriptions, and sorts the resulting data frame.
}
\examples{
# Assuming `model` is a preGroup model
\dontrun{
set.seed(42)
mod <- preGroup(Species ~ ., data = iris, treatment_indicator = "Petal.Width")
result <- coef(mod)
result
}
}
