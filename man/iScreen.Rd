\name{iScreen}
\alias{iScreen}
\title{image-based HTS analysis}
\usage{
  iScreen(data = NULL, formula = NULL, control = NULL,
    FUN = NULL, ...)
}
\arguments{
  \item{data}{A data frame containing the variables in the
  model with the first column as Well ID. For details of
  the Well ID format, see \code{\link{WellCheck}}. If not
  found, an error will be reported.}

  \item{formula}{An object of class "\code{\link{formula}}"
  (or one that can be coerced to that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under 'Details' from
  \code{\link{glm}}.}

  \item{control}{A optional vector of the length as the
  number of rows in data indicating which row or rows from
  data should be used as control when fitting generalized
  linear model.}

  \item{FUN}{A user-provided function for analysis}

  \item{...}{Arguments to be passed to method. For details
  also see \code{\link{formula}}.}
}
\value{
  A "iScreen" list of two components, fit and coefficients.
  \item{fit}{A fitted generalized linear model object from
  \code{\link{glm}} or NULL if FUN is used}
  \item{coefficients}{Estimations of coefficients for each
  Well ID}
}
\description{
  Analysis of image-based high-throughput RNAi screen via
  either generalized linear model or customized user
  function.
}
\details{
  For argument FUN, user has to provide a custom function.
  Function argument has to be from data. Return value can
  be a single value or a vector. For vector return, the
  first value has to be estimate for each Well ID.
}
\examples{
data(autophagy)
fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control  == 1))
head(fit.auto$coefficients)
}

