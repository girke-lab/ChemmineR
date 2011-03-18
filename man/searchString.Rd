\name{searchString}
\alias{searchString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PubChem Similarity (Fingerprint) SMILES Search 
}
\description{
Accepts one SMILES string (Simplified Molecular Input Line Entry Specification)
and performs a >0.95 similarity PubChem fingerprint search, returning the 
hits in an \code{SDFset} container. The ChemMine Tools web service
is used as an intermediate, to translate queries from plain HTTP POST to
a PubChem Power User Gateway (PUG) query.
}
\usage{
searchString(smiles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smiles}{
  	  A \code{character} object which contains one SMILES string
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
\dontrun{
## search a compound on PubChem
compounds <- searchString("CC(=O)OC1=CC=CC=C1C(=O)O")}
}

\keyword{ utilities }
