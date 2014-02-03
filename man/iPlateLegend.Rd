\name{iPlateLegend}
\alias{iPlateLegend}
\title{Plotting iPlate legend}
\usage{
  iPlateLegend(object, ...)
}
\arguments{
  \item{object}{A iPlate object.}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{image}} and \code{\link{par}}.}
}
\description{
  Function for plotting legend for object returned by
  iPlate. For more details about the graphical parameter
  arguments, see \code{\link{image}} and \code{\link{par}}.
}
\examples{
data(autophagy)
p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
iPlateLegend(p1)
}

