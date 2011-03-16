\name{getIds}
\alias{getIds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
numeric to \code{SDFset}
}
\description{
Accepts one or more PubChem compound ids
and downloads the corresponding compounds from PubChem Power User Gateway (PUG)
returning results in an \code{SDFset} container. The ChemMine Tools web service
is used as an intermediate, to translate queries from plain HTTP POST to
a PUG SOAP query.
}
\usage{
getIds(cids)
}

\arguments{
  \item{cids}{
  	  A \code{numeric} object which contains one or more PubChem cids
}
}
\value{
\item{SDFset}{ for details see ?"SDFset-class"}
}
\references{
PubChem PUG SOAP: http://pubchem.ncbi.nlm.nih.gov/pug_soap/pug_soap_help.html

Chemmine web service: http://chemmine.ucr.edu

PubChem help: http://pubchem.ncbi.nlm.nih.gov/search/help_search.html
}
\author{
Tyler Backman
}

\examples{
## fetch 2 compounds from PubChem
compounds <- getIds(c(111,123))
}

\keyword{ utilities }
