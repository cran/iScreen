\name{iScreenPlot}
\alias{iScreenPlot}
\title{Plotting iScreen}
\usage{
  iScreenPlot(object, xlab = "Column", ylab = "Row", ...)
}
\arguments{
  \item{object}{A iScreen object.}

  \item{xlab}{Default is "Column".}

  \item{ylab}{Default is "Row".}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{plot}} and \code{\link{par}}}
}
\description{
  Function for plotting object returned by iScreen. For
  more details about the graphical parameter arguments, see
  \code{\link{par}}.
}
\examples{
data(autophagy)
fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
iScreenPlot(fit.auto)
}

