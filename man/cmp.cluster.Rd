\name{cmp.cluster}
\alias{cmp.cluster}
\title{cluster compounds using a descriptor database}
\description{
	'cmp.cluster' uses structural compound descriptors and clusters the
		compounds based on their pairwise distances. \code{cmp.cluster} uses
		single linkage to measure distance between clusters when it
		merges clusters. It accepts both a single cutoff and a
		cutoff vector. By using a cutoff vector, it can generate results
		similar to hierarchical clustering after tree cutting.
}
\usage{
cmp.cluster(db, cutoff, is.similarity = TRUE, save.distances = FALSE,
        use.distances = NULL, quiet = FALSE, ...)
}
\arguments{
  \item{db}{The desciptor database, in the format returned by 'cmp.parse'.}
  \item{cutoff}{The clustering cutoff. Can be a single value or a vector. The
      cutoff gives the maximum distance between two compounds in order to
      group them in the same cluster.}
  \item{is.similarity}{Set when the cutoff supplied is a similarity cutoff.
      This cutoff is the minimum similarity value between two compounds such
      that they will be grouped in the same cluster.}
  \item{save.distances}{whether to save distance for future clustering. See
	  details below.}
  \item{use.distances}{Supply pre-computed distance matrix.}
  \item{quiet}{Whether to suppress the progress information.}
  \item{...}{Further arguments to be passed to \code{cmp.similarity}.}
}
\details{
    \code{cmp.cluster} will compute distances on the fly if \code{use.distances} is not set.
    Furthermore, if \code{save.distances} is not set, the distance values computed will never be
    stored and any distance between two compounds is guaranteed not to be
    computed twice. Using this method, \code{cmp.cluster} can deal with large databases
    when a distance matrix in memory is not feasible. The speed of the clustering
    function should be slowed when using a transient distance calculation.
	
    When \code{save.distances} is set, \code{cmp.cluster} will be forced to compute the
    distance matrix and save it in memory before the clustering. This is
    useful when additional clusterings are required in the future without re-computed 
    the distance matrix. Set \code{save.distances} to TRUE if you
    only want to force the clustering to use this 2-step approach; otherwise,
    set it to the filename under which you want the distance matrix to be
    saved. After you save it, when you need to reuse the distance matrix, you
    can 'load' it, and supply it to \code{cmp.cluster} via the \code{use.distances} argument.

    \code{cmp.cluster} supports a vector of several cutoffs. When you have multiple cutoffs,
    \code{cmp.cluster} still guarantees that pairwise distances will never be
    recomputed, and no copy of distances is kept in memory. It is guaranteed to
    be as fast as calling \code{cmp.cluster} with a single cutoff that results in the
    longest processing time, plus some small overhead linear in processing
    time.
}
\value{
    Returns a \code{data.frame}. Besides a variable giving compound ID, each of the
    other variables in the data frame will either give the cluster IDs of
    compounds under some clustering cutoff, or the size of clusters that the
    compounds belong to. When N cutoffs are given, in total 2*N+1 variables
    will be generated, with N of them giving the cluster ID of each compound
    under each of the N cutoffs, and the other N of them giving the cluster
    size under each of the N cutoffs. The rows are sorted by cluster sizes.
}
\author{Y. Eddie Cao, Li-Chang Cheng}
\seealso{\code{\link{cmp.parse1}}, \code{\link{cmp.parse}}, \code{\link{cmp.search}}, \code{\link{cmp.similarity}}}
\examples{
## Load sample SD file
# data(sdfsample); sdfset <- sdfsample

## Generate atom pair descriptor database for searching
# apset <- sdf2ap(sdfset) 

## Loads atom pair and atom pair fingerprint samples provided by library
data(apset) 
db <- apset
fpset <- desc2fp(apset)

## Clustering of 'APset' object with multiple cutoffs
clusters <- cmp.cluster(db=apset, cutoff=c(0.5, 0.85))

## Clustering of 'FPset' object with multiple cutoffs. This method allows to call 
## various similarity methods provided by the fpSim function.
clusters2 <- cmp.cluster(fpset, cutoff=c(0.5, 0.7), method="Tversky") 

## Saves the distance matrix before clustering:
clusters <- cmp.cluster(db, cutoff=0.65, save.distances="distmat.rda")
# Later one reload the matrix and pass it the clustering function. 
load("distmat.rda")
clusters <- cmp.cluster(db, cutoff=0.60, use.distances=distmat)
}
\keyword{utilities}
