\name{searchSub}
\alias{searchSub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{SDFset} to \code{SDFset}
}
\description{
Accepts one \code{SDFset} container
and performs a >0.95 similarity PubChem fingerprint search, returning the 
hits in an \code{SDFset} container. The ChemMine Tools web service
is used as an intermediate, to translate queries from plain HTTP POST to
a PubChem Power User Gateway (PUG) query. If the input object
contains multiple items, only the first is used as a query.
}
\usage{
searchSub(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
  	  A \code{SDFset} object which contains one compound
}
}
\value{
\item{SDFset}{ for details see ?"SDFset-class"}
}
\references{
PubChem PUG SOAP: http://pubchem.ncbi.nlm.nih.gov/pug_soap/pug_soap_help.html

Chemmine web service: http://chemmine.ucr.edu

PubChem help: http://pubchem.ncbi.nlm.nih.gov/search/help_search.html

SMILES Format: http://en.wikipedia.org/wiki/Chemical_file_format#SMILES
}
\author{
Tyler Backman
}

\examples{
## get a sample compound
data(sdfsample); sdfset <- sdfsample[1]
## search a compound on PubChem
compounds <- searchSub(sdfset)
}

\keyword{ utilities }
