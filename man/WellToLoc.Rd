\name{WellToLoc}
\alias{WellToLoc}
\title{Convert WellID}
\usage{
  WellToLoc(WellNO)
}
\arguments{
  \item{WellNO}{A character or a vector of characters as
  variable WellID to converted}
}
\value{
  A list consisting of 4 vectors. \item{x.loc}{A vector of
  the length of WellNO indicating x coordinate of each
  WellNO upon plot} \item{y.loc}{A vector of the length of
  WellNO indicating y coordinate of each WellNO upon plot}
  \item{x label}{A vector of names to be used as label on X
  axis upon plot} \item{y.label}{A vector of names to be
  used as label on Y axis upon plot}
}
\description{
  Convert Well ID into plot coordinate and label
}

