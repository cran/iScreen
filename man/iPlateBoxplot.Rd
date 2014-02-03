\name{iPlateBoxplot}
\alias{iPlateBoxplot}
\title{Plotting iPlate boxplot}
\usage{
  iPlateBoxplot(object, by = "row", ...)
}
\arguments{
  \item{object}{A iPlate object.}

  \item{by}{Default is "row" and can also be "column".}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{boxplot}} and \code{\link{par}}.}
}
\description{
  Function for box plotting object returned by iPlate. For
  more details about the graphical parameter arguments, see
  \code{\link{boxplot}} and \code{\link{par}}.
}
\examples{
data(autophagy)
p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
iPlateBoxplot(p1)
}

