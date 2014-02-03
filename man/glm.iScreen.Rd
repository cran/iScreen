\name{glm.iScreen}
\alias{glm.iScreen}
\title{Generalized linear model for image-based HTS}
\usage{
  glm.iScreen(data, formula, control = NULL, ...)
}
\arguments{
  \item{data}{A data frame containing the variables in the
  model with the first column as Well ID. For details of
  the Well ID format, see \code{\link{WellCheck}}. If not
  found, an error will be reported.}

  \item{formula}{An object of class \code{\link{formula}}or
  one that can be coerced to that class. a symbolic
  description of the model to be fitted. The details of
  model specification are given under 'Details' from
  \code{\link{glm}}.}

  \item{control}{A optional vector of the length as the
  number of rows in data indicating which row or rows from
  data should be used as control when fitting generalized
  linear model.}

  \item{...}{Arguments to be passed to method. For details
  also see \code{\link{formula}}.}
}
\value{
  A list of two compoments, fit and coefficients.
  \item{fit}{A fitted generlized linear model object from
  \code{\link{glm}}} \item{coefficients}{Estimations of
  coefficients for each Well ID}
}
\description{
  Perform generalized linear model for image-based HTS data
  analysis.
}

