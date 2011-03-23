\name{atomprop}
\alias{atomprop}
\docType{data}
\title{
Standard atomic weights
}
\description{
Data frame with atom names, symbols, standard atomic weights, group number and period number.
}
\usage{data(atomprop)}
\format{
  The format is a data frame with 117 rows and 6 columns.
}
\source{
Columns 1 to 4 from: http://iupac.org/publications/pac/78/11/2051/
Columns 5 to 6 from: http://en.wikipedia.org/wiki/List_of_elements
}
\references{
Pure Appl. Chem., 2006, Vol. 78, No. 11, pp. 2051-2066
}
\examples{
data(atomprop)
atomprop[1:4,]
}
\keyword{datasets}
