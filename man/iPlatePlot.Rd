\name{iPlatePlot}
\alias{iPlatePlot}
\title{Plotting iPlate object}
\usage{
  iPlatePlot(object, ...)
}
\arguments{
  \item{object}{A iPlate object.}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{plot}} and \code{\link{par}}}
}
\description{
  Function for plotting object returned by iPlate. For more
  details about the graphical parameter arguments, see
  \code{\link{par}}.
}
\examples{
data(autophagy)
p1 <- iPlate(autophagy, "dot.number", log=TRUE, FUN=mean)
iPlatePlot(p1)
}

