\name{WellCheck}
\alias{WellCheck}
\title{Check WellID format}
\usage{
  WellCheck(x)
}
\arguments{
  \item{x}{A character or a vector of characters as
  variable WellID}
}
\value{
  The value is logical of length one. Return TRUE if and
  only if all elements in x accord with WELLID requirement,
  otherwise return FALSE.
}
\description{
  Examination of variable WellID format
}
\examples{
#WellID format should be "A01", "B02", "C03", ..."J12", ...
WellCheck(c("A01"))
WellCheck(c("A1"))
WellCheck(c("A01", "B02", "C03"))
WellCheck(c("A01", "B2", "C03"))
}

