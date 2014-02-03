\name{iScreenLegend}
\alias{iScreenLegend}
\title{Plotting iScreen legend}
\usage{
  iScreenLegend(object, ...)
}
\arguments{
  \item{object}{A iScreen object.}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{image}} and \code{\link{par}}.}
}
\description{
  Function for plotting legend for object returned by
  iScreen. For more details about the graphical parameter
  arguments, see \code{\link{image}} and \code{\link{par}}.
}
\examples{
data(autophagy)
fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
iScreenLegend(fit.auto)
}

