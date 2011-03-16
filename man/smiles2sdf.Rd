\name{smiles2sdf}
\alias{smiles2sdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{character} to \code{SDFset}
}
\description{
Accepts one compound as a
SMILES string (Simplified Molecular Input Line Entry Specification) and
returns it's equivalent as an \code{SDFset} container. 
The compound is submitted to the ChemMine Tools web service for conversion
with the Open Babel Open Source Chemistry Toolbox.
}
\usage{
smiles2sdf(smiles)
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
Chemmine web service: http://chemmine.ucr.edu

Open Babel: http://openbabel.org

SMILES Format: http://en.wikipedia.org/wiki/Chemical_file_format#SMILES
}
\author{
Tyler Backman
}

\examples{
## convert to sdf
sdf <- smiles2sdf("CC(=O)OC1=CC=CC=C1C(=O)O\tname")
}
\keyword{ utilities }
