\name{jarvisPatrick}
\alias{jarvisPatrick}
\title{
Jarvis-Patrick Clustering
}
\description{
Function to perform Jarvis-Patrick clustering. The algorithm requires a nearest neighbor table, which consists of \code{j} nearest neighbors for each item in the dataset. This information is then used to join items into clusters that share at least \code{k} nearest neighbors. The values for \code{j} and \code{k} are user-defined parameters. The \code{jarvisPatrick} function can generate the nearest neighbor table for \code{APset} and \code{FPset} objects and then perform Jarvis-Patrick clustering on that table. It also accepts a precomputed nearest neighbor table in form of an object of class matrix. The output is a cluster vector with the item labels in the name slot and the cluster IDs in the data slot. Alternatively, the function can return the nearest neighbor matrix. As third parameter the user can set a minimum similarity \code{cutoff} value for generating the nearest neighbor table. The latter is an optional setting that is not part of the original Jarvis-Patrick algorithm. It allows to generate more tight clusters and to minimize some limitations of this method, such as joining completely unrelated items when clustering small data sets.
}
\usage{
jarvisPatrick(x, j, k, cutoff = NA, type = "cluster", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Input of class \code{FPset}, \code{APset} or \code{matrix}. 
}
  \item{j}{
Number of nearest neighbors to consider when generating nearest neighbor table.
}
  \item{k}{
Minimum number of nearest neighbors two rows (items) in nearest neighbor table need to have in common to join them into the same cluster. 
}
  \item{cutoff}{
Uses a minimum similarity threshold when generating the nearest neighbor table. This parameter allows to obtain tighter clustering results. It is optional and not part of the original Jarvis-Patrick algorithm. Under the default settings this behavior is turned off. If \code{cutoff} is set to a stringent value then it is often required to set \code{k} to a small value.  
}
  \item{type}{
If \code{type = "cluster"}, Jarvis-Patrick clustering is performed and the result is returned in a named vector. If \code{type = "matrix"}, the nearest neighbors table is returned as \code{matrix}. 
}
  \item{\dots}{
Arguments to be passed to/from other methods.
}
}
\details{
...
}
\value{
Depending on the setting under the \code{type} argument, the function returns the clustering result in a \code{named vector} or a nearest neighbor table as \code{matrix}.
}
\references{
Jarvis RA, Patrick EA (1973) Clustering Using a Similarity Measure Based on Shared Near Neighbors. IEEE Transactions on Computers, C22, 1025-1034. URLs: http://davide.eynard.it/teaching/2012_PAMI/JP.pdf, http://www.daylight.com/dayhtml/doc/cluster/index.pdf
}
\author{
Thomas Girke
}
\note{
...
}
\seealso{
Functions: \code{cmp.cluster}
}
\examples{
## Load/create sample APset and FPset 
data(apset)
fpset <- desc2fp(apset)

## Standard Jarvis-Patrick clustering on APset/FPset objects
jarvisPatrick(apset, j=6, k=5)
jarvisPatrick(fpset, j=6, k=5)

## Modified Jarvis-Patrick clustering with minimum similarity 'cutoff' 
## value (here Tanimoto coefficient)
jarvisPatrick(fpset, j=6, k=2, cutoff=0.6, method="Tanimoto")

## Output nearest neighbor table (matrix)
nnm <- jarvisPatrick(fpset, j=6, k=5, type="matrix")

## Perform clustering on precomputed nearest neighbor table
jarvisPatrick(nnm, j=6, k=5)
}
\keyword{utilities}
