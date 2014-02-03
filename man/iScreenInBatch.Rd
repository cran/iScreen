\name{iScreenInBatch}
\alias{iScreenInBatch}
\title{Batch process of image-based HTS}
\usage{
  iScreenInBatch(data = NULL, formula = NULL,
    control = NULL, FUN = NULL, ...)
}
\arguments{
  \item{data}{A data frame containing the variables in the
  model with the first two column as Well ID and Plate ID.
  For details of the Well ID format, see
  \code{\link{WellCheck}}. If not found, an error will be
  reported.}

  \item{formula}{An object of class "\code{\link{formula}}"
  (or one that can be coerced to that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under 'Details' from
  \code{\link{glm}}.}

  \item{control}{A optional verctor of the length as the
  number of rows in data indicating which row or rows from
  data should be used as control when fitting generalized
  linear model.}

  \item{FUN}{A user-provided function for anlaysis.}

  \item{...}{Arguments to be passed to method. For details
  also see \code{\link{formula}}.}
}
\value{
  A "iScreenInBatch" data frame containing summary
  coefficients.
}
\description{
  Batch processing of image-based high-throughput RNAi
  screen via either generalized linear model or customized
  user function.
}
\details{
  For argument FUN, user has to providde a custom function.
  Function argument has to be from data. Return value can
  be a single value or a vector. For vector return, the
  first value has to be estimate for each Well ID.
}

