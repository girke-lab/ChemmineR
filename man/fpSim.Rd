\name{fpSim}
\alias{fpSim}
\title{
Fingerprint Search 
}
\description{
Search function for fingerprints, such as PubChem or atom pair fingerprints. Enables structure similarity comparisons, searching and clustering.
}
\usage{
fpSim(x, y, sorted=TRUE, method="Tanimoto", cutoff=0, top="all", alpha=1, beta=1, ...)
}
\arguments{
  \item{x}{
Query molecule of class \code{numeric}, \code{FP} or \code{FPset} (of length one) containing binary fingerprint data. Both \code{x} and \code{y} need to have the same number of bits and should contain the same type of fingerprints.
}
  \item{y}{
Subject molecule(s) of class \code{numeric}, \code{matrix}, \code{FP} or \code{FPset} containing binary fingerprint data.
}
  \item{sorted}{
return results sorted or unsorted
}
  \item{method}{
Similarity coefficient to return. One can choose here from several predefined similarity measures: "Tanimoto" (default), "Euclidean", "Tversky" or "Dice". Alternatively, one can pass on any custom similarity function containing the arguments a, b, c and d. For instance, one can define "myfct <- function(a, b, c, d) c/(alpha*a + beta*b + c)" and then pass on \code{method=myfct}. The variable 'c' is the number of "on-bits" common in both compounds, 'd' is the number of "off-bits" common in both compounds, and 'a' and 'b' are the number of "on-bits" that are unique in one or the other compound, respectively.
}
  \item{cutoff}{
allows to restrict results to hits above a similarity cutoff value; default \code{cutoff=0} returns results for all compounds in \code{y}.
}
  \item{top}{
allows to restrict number of subject molecules to return; default \code{top="all"} returns results for all compounds in \code{y} above \code{cutoff} value.
}
  \item{alpha}{
Only used when method="Tversky". Allows to specify the weighting variable 'alpha' of the Tversky index: c/(alpha*a + beta*b + c)
}
  \item{beta}{
Only used when method="Tversky". Allows to specify the weighting variable 'beta' of the Tversky index. 
}
\item{\dots}{
arguments to be passed to/from other methods.
}
}
\details{
...
}
\value{
Returns \code{numeric vector} with similarity coefficients as values and compound identifiers as names. 
}
\references{
Tanimoto similarity coefficient: Tanimoto TT (1957) IBM Internal Report 17th Nov see also Jaccard P (1901) Bulletin del la Societe Vaudoisedes Sciences Naturelles 37, 241-272.

PubChem fingerprint specification: ftp://ftp.ncbi.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
}
\author{
Thomas Girke
}
\note{
...
}

\seealso{
Functions: \code{fp2bit} 
}
\examples{
## Load PubChem SDFset sample
data(sdfsample); sdfset <- sdfsample
cid(sdfset) <- sdfid(sdfset)

## Convert base 64 encoded fingerprints to character vector or binary matrix
fpset <- fp2bit(sdfset)

## Alternatively, one can use atom pair fingerprints 
\dontrun{
fpset <- desc2fp(sdf2ap(sdfset))
}

## Pairwise compound structure comparisons
fpSim(x=fpset[1], y=fpset[2], method="Tanimoto")

## Structure similarity searching: x is query and y is fingerprint database  
fpSim(x=fpset[1], y=fpset) 

## Controlling the output
fpSim(x=fpset[1], y=fpset, method="Tversky", cutoff=0.4, top=4, alpha=0.5, beta=1) 

## Use custom distance function
myfct <- function(a, b, c, d) c/(a+b+c+d)
fpSim(x=fpset[1], y=fpset, method=myfct) 

## Compute fingerprint-based Tanimoto similarity matrix 
simMA <- sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) 

## Hierarchical clustering with simMA as input
hc <- hclust(as.dist(1-simMA), method="single")

## Plot hierarchical clustering tree
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE) 
}
\keyword{utilities}
