\name{iWell}
\alias{iWell}
\title{image-based HTS Well object}
\usage{
  iWell(x, y, d, c, n = 4, angle = 0, type = 1)
}
\arguments{
  \item{x}{A vector containing X coordinate of each plot
  unit.}

  \item{y}{A vector containing Y coordinate of each plot
  unit.}

  \item{d}{A vector containing diameter of each plot unit.}

  \item{c}{A vector containing plot colore of each plot
  unit.}

  \item{n}{A numeric or vector containing number of sides
  for each plot unit.}

  \item{angle}{A numeric or vector containing rotation
  angle of each plot unit, in degrees.}

  \item{type}{A numeric or vector containing plotting type
  type=1 => interior filled, type=2 => edge, type=3 =>
  both.}
}
\value{
  A object "iWell".
}
\description{
  image-based HTS Well object to be used for visualization
}
\details{
  For details also see \code{\link{ngon}}.
}

