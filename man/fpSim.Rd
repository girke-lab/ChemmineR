\name{fpSim}
\alias{fpSim}
\title{
PubChem Fingerprint Search 
}
\description{
Function to use PubChem fingerprints for structure similarity comparisons, searching and clustering.
}
\usage{
fpSim(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{vector} containing binary fingerprint data. Needs to have the same length as \code{y} (\code{vector} or \code{matrix} row).
}
  \item{y}{
\code{vector} or \code{matrix} containing binary fingerprint data.
}
}
\details{
The function computes the Tanimoto coefficients for pairwise comparisons of binary fingerprints. The coefficient is defined as c/(a+b+c), which is the proportion of the "on-bits" shared among the fingerprints of two compounds divided by their union. The variable c is the number of "on-bits" common in both compounds, while a and b are the number of "on-bits" that are unique in one or the other compound, respectively.

}
\value{
Returns \code{numeric vector} with Tanimoto coefficients as values and compound identifiers as names. 
}
\references{
Tanimoto similarity coefficient: Tanimoto TT (1957) IBM Internal Report 17th Nov see also Jaccard P (1901) Bulletin del la Societe Vaudoisedes Sciences Naturelles 37, 241-272.

PubChem fingerprint specification: ftp://ftp.ncbi.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
}
\author{
Thomas Girke
}
\note{
Limitation: PubChem fingerprints need to be provided, such as in PubChem's SD files.  
}

\seealso{
Functions: \code{fp2bit} 
}
\examples{
## Load PubChem SDFset sample
data(sdfsample); sdfset <- sdfsample
cid(sdfset) <- sdfid(sdfset)

## Convert base 64 encoded fingerprints to character vector or binary matrix
fpset <- fp2bit(x=sdfset, type=1)
fpset <- fp2bit(x=sdfset, type=2)

## Pairwise compound structure comparisons
fpSim(x=fpset[1,], y=fpset[2,]) 

## Structure similarity searching: x is query and y is fingerprint database  
fpSim(x=fpset[1,], y=fpset) 

## Compute fingerprint-based Tanimoto similarity matrix 
simMA <- sapply(rownames(fpset), function(x) fpSim(x=fpset[x,], fpset)) 

## Hierarchical clustering with simMA as input
hc <- hclust(as.dist(simMA), method="single")

## Plot hierarchical clustering tree
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE) 
}
\keyword{utilities}
