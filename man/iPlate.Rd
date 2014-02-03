\name{iPlate}
\alias{iPlate}
\title{image-based HTS Plate}
\usage{
  iPlate(data, column, log = FALSE, FUN = mean, ...)
}
\arguments{
  \item{data}{A data frame containing the variables in the
  model with the first column as Well ID. For details of
  the Well ID format, see \code{\link{WellCheck}}. If not
  found, an error will be reported.}

  \item{column}{Specify the column to be used for analysis
  and visualization.}

  \item{log}{Specify if data need be logarithm
  transformed.}

  \item{FUN}{A function to compute the summary statistics
  which can be applied to all data subsets.}

  \item{...}{Further arguments passed to or used by
  methods. See \code{\link{aggregate}} for further
  information.}
}
\value{
  A object "iPlate" list of 4 elements. \item{z}{A data
  frame containing Well ID and summary statistics for each
  Well ID} \item{loc}{A list of coordinate information
  returned by \code{\link{WellToLoc}}} \item{log}{A logical
  value indicating if data is logarithm transformed}
  \item{FUN}{the function to compute the summary
  statistics}
}
\description{
  image-based HTS Plate object to be used for analysis and
  visualization
}

