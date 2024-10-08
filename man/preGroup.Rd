% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preGroup.r
\name{preGroup}
\alias{preGroup}
\title{Pre-Group Function for Building Ensemble Models}
\usage{
preGroup(
  formula,
  data,
  treatment_indicator,
  alpha.mvs = c(0, 1),
  family = "gaussian",
  ad.alpha = NA,
  ad.penalty = "lambda.min",
  use.grad = TRUE,
  weights,
  type = "both",
  sampfrac = 0.5,
  maxdepth = 3L,
  learnrate = 0.01,
  mtry = Inf,
  ntrees = 500,
  confirmatory = NULL,
  singleconditions = FALSE,
  winsfrac = 0.025,
  normalize = TRUE,
  standardize = FALSE,
  ordinal = TRUE,
  nfolds = 10L,
  tree.control,
  tree.unbiased = TRUE,
  removecomplements = TRUE,
  removeduplicates = TRUE,
  verbose = FALSE,
  par.init = FALSE,
  par.final = FALSE,
  sparse = FALSE,
  ...
)
}
\arguments{
\item{formula}{an object of class 'formula' describing the model to be fitted.}

\item{data}{a data frame containing the variables specified in the formula.}

\item{treatment_indicator}{a character string specifying the column used for treatment division.}

\item{alpha.mvs}{a dummay vector of length two, indicating the alpha values used
within  and between group penalization. A value of 1 yields lasso, a value of 0 ridge, a value
between 0 and 1 elastic net. The default uses ridge for feature selection and lasso
for group selection.}

\item{family}{a description of the error distribution to be used in the model. It can be one
of 'gaussian', 'binomial', etc.}

\item{ad.alpha}{optional alpha parameter for adjustment.}

\item{ad.penalty}{the penalty type to be used in adjustment, default is 'lambda.min'.}

\item{use.grad}{logical, indicating whether to use gradients.}

\item{weights}{an optional vector of 'weights' to be used in the fitting process.}

\item{type}{the type of model components to include: 'rules', 'linear', or 'both'.}

\item{sampfrac}{the fraction of the data to be used in each bootstrap sample.}

\item{maxdepth}{maximum depth of the trees in the model.}

\item{learnrate}{learning rate for boosting.}

\item{mtry}{number of variables randomly sampled as candidates at each split.}

\item{ntrees}{number of trees to grow.}

\item{confirmatory}{optional vector of confirmatory conditions.}

\item{singleconditions}{logical, to include only single conditions in rules.}

\item{winsfrac}{the fraction of extreme values to winsorize.}

\item{normalize}{logical, whether to normalize predictor variables.}

\item{standardize}{logical, whether to standardize predictor variables.}

\item{ordinal}{logical, to treat ordered factors as numeric.}

\item{nfolds}{number of folds for cross-validation.}

\item{tree.control}{list of control parameters for tree construction.}

\item{tree.unbiased}{logical, indicating whether to use unbiased trees (conditional inference trees; ctrees) or biased trees (recursive partitionaing; rpart).}

\item{removecomplements}{logical, whether to remove complement rules.}

\item{removeduplicates}{logical, whether to remove duplicate rules.}

\item{verbose}{logical, indicating whether to print detailed output during fitting.}

\item{par.init}{logical, indicating whether to parallelize the initial fit.}

\item{par.final}{logical, indicating whether to parallelize the final model fitting.}

\item{sparse}{logical, to use sparse model matrices.}

\item{...}{additional arguments affecting the model fitting.}
}
\value{
A list containing the fitted model object, the original call, and a classification
of the model rules into types such as 'linear', 'prognostic', and 'prescriptive'.
}
\description{
This function facilitates the creation of ensemble models by grouping data based on
treatment indicators and other factors. It allows for various configurations and customizations
to tailor the model fitting process to specific needs, including handling of different
statistical family distributions.
}
\examples{
\dontrun{
 data(iris)
 result <- preGroup(Species ~ ., data = iris, treatment_indicator = "Petal.Width")
 result
}

}
