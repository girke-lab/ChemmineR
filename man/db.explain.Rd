\name{db.explain}
\alias{db.explain}
\title{Explain an atom-pair descriptor or an array of atom-pair descriptors}
\description{
'db.explain' will take an atom-pair descriptor in numeric or a set of such descriptors, and interpret what they represent in a more human readable way.
}
\usage{
db.explain(desc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{desc}{The descriptor or the array/vector of descriptors}
}
\details{
	'desc' can be a single numeric giving a single descriptor or can be
	any container data type, such as vector or array, such that
	'length(desc)' returns 2 or larger.
}
\value{
    Return a character vector describing the descriptors.	
}

\seealso{\code{\link{cmp.parse}}}
\examples{
## Load sample SD file
# data(sdfsample); sdfset <- sdfsample

## Generate atom pair descriptor database for searching
# apset <- sdf2ap(sdfset) 

## Loads same atom pair sample data set provided by library
data(apset) 
db <- apset

## Return atom pairs of first compound in human readable format
db.explain(db[1])
}
\keyword{utilities}
