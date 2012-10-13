\name{fp2bit}
\alias{fp2bit}
\title{
Convert base 64 fingerprints to binary
}
\description{
The function converts the base 64 encoded PubChem fingerprints to a binary \code{matrix} or a \code{character} vector. If applied to a \code{SDFset} object, then its data block needs to contain the PubChem fingerprint information.
}
\usage{
fp2bit(x, type = 3, fptag = "PUBCHEM_CACTVS_SUBSKEYS")
}
\arguments{
  \item{x}{
Object of class \code{SDFset}, \code{matrix} or \code{character}
}
  \item{type}{
If set to \code{1}, the results are returned as binary \code{matrix}. If set to \code{2}, the results are returned as \code{character} strings in a named \emph{vector}. If set to \code{3} (default), the results are returned as \code{FPset} object. 
}
  \item{fptag}{
Name tag in SDF data block where the PubChem fingerprints are stored. Default is set to "PUBCHEM_CACTVS_SUBSKEYS".
}
}
\details{
...
}
\value{
\code{matrix}, \code{character} or \code{FPset}
}
\references{
See PubChem fingerprint specification at: ftp://ftp.ncbi.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
}
\author{
Thomas Girke
}
\seealso{
Functions: \code{fpSim}
}
\examples{
## Load PubChem SDFset sample
data(sdfsample); sdfset <- sdfsample
cid(sdfset) <- sdfid(sdfset)

## Convert base 64 encoded fingerprints to FPset object
fpset <- fp2bit(sdfset)

## Pairwise compound structure comparisons
fpSim(fpset[1], fpset[2]) 

## Structure similarity searching: x is query and y is fingerprint database
fpSim(x=fpset[1], y=fpset, method="Tanimoto", cutoff=0, top="all") 

## Compute fingerprint based Tanimoto similarity matrix 
simMA <- sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE)) 

## Hierarchical clustering with simMA as input
hc <- hclust(as.dist(1-simMA), method="single")

## Plot hierarchical clustering tree
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE)
}
\keyword{utilities}
