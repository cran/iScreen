\name{iScreenBoxplot}
\alias{iScreenBoxplot}
\title{Plotting iScreen boxplot}
\usage{
  iScreenBoxplot(object, by = "row", ...)
}
\arguments{
  \item{object}{A iScreen object.}

  \item{by}{Default is "row" and can also be "column".}

  \item{...}{Arguments to be passed to methods. See
  \code{\link{boxplot}} and \code{\link{par}}.}
}
\description{
  Function for box plotting object returned by iScreen. For
  more details about the graphical parameter arguments, see
  \code{\link{boxplot}} and \code{\link{par}}.
}
\examples{
data(autophagy)
fit.auto <- iScreen(autophagy, dot.number~WellID, family=poisson, control=(autophagy$control == 1))
iScreenBoxplot(fit.auto)
}

