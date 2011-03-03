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
# load sample database from web
db <- cmp.parse("http://bioweb.ucr.edu/ChemMineV2/static/example_db.sdf")
# explain descriptor 1 of compound 1
db.explain(db$descdb[[1]][[1]])
# explain descriptor 1 to 10 of compound 1
db.explain(db$descdb[[1]][1:10])
# explain all descriptors of compound 1
db.explain(db$descdb[[1]])
}
\keyword{utilities}
