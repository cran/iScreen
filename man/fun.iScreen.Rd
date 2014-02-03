\name{fun.iScreen}
\alias{fun.iScreen}
\title{Custom analysis of image-based HTS}
\usage{
  fun.iScreen(data, FUN)
}
\arguments{
  \item{data}{A data frame containing the variables in the
  model with the first column as Well ID. For details of
  the Well ID format, see \code{\link{WellCheck}}. If not
  found, an error will be reported.}

  \item{FUN}{A user-provided function for analysis}
}
\value{
  A list of two components, fit and coefficients.
  \item{fit}{NULL} \item{coefficients}{Estimations of
  coefficients for each Well ID}
}
\description{
  Perform custom analysis of image-based HTS data.
}
\details{
  For argument FUN, user has to provide a custom function.
  Function argument has to be from data. Return value can
  be a single value or a vector. For vector return, the
  first value has to be estimate for each Well ID.
}

