\name{cmp.cluster}
\alias{cmp.cluster}
\title{cluster compounds using a descriptor database}
\description{
	'cmp.cluster' uses compound descriptors in a database and clusters these
		compounds based on their pairwise distances. 'cmp.cluster' uses
		single linkage to measure distance between clusters when it
		merges clusters. 'cmp.cluster' accepts both a single cutoff and a
		cutoff vector. By using a cutoff vector, it can generate the
		same result as hierachical clustering.
}
\usage{
cmp.cluster(db, cutoff, is.similarity = TRUE, save.distances = FALSE,
        use.distances = NULL, quiet = FALSE, ...)
}
\arguments{
  \item{db}{The desciptor database, in the format returned by 'cmp.parse'.}
  \item{cutoff}{The clustering cutoff. Can be a single value or a vector. The
      cutoff gives the maximum distance between two compounds in order to
      group them in the same clsuter.}
  \item{is.similarity}{Set when the cutoff supplied is a similarity cutoff.
      This cutoff is the mimumum similarity value between two compounds such
      that they will be grouped in the same cluster.}
  \item{save.distances}{whether to save distance for future clustering. See
	  details below.}
  \item{use.distances}{Supply pre-computed distance matrix.}
  \item{quiet}{Whether to supress the progress information.}
  \item{...}{Further arguments to be passed to 'cmp.similarity' to calculate similarities if necessary.}
}
\details{
    'cmp.cluster' will compute distance on the fly if 'use.distances' is not set.
    Furthermore, if 'save.distances' is not set, the distance will never be
    stored and distance between any two compounds is guaranteed not to be
    computed twice. Using this method, 'cmp.cluster' can deal with large database,
    when a distance matrix in memory is not feasible. The speed of this cluster
    function should be slowed because of using this transient distance value.
	
    When 'save.distances' is set, 'cmp.cluster' will be forced to compute the
    distance matrix and save it in memory before doing clustering. This is
    useful when you need to do further clustering in the future and do not want
    the distance to be re-computed then. Set 'save.distances' to TRUE if you
    only want to force the clustering to use this 2-step approach; otherwise,
    set it to the filename under which you want the distance matrix to be
    saved. After you save it, when you need to reuse the distance matrix, you
    can 'load' it, and supply to 'cmp.cluster' via the 'use.distances' argument.

    'cmp.cluster' supports vector of cutoffs. When you have multiple cutoffs,
    'cmp.cluster' still guarantees that pairwise distances will never be
    recomputed, and no copy of distances is kept in memory. It is guaranteed to
    be as fast as calling 'cmp.cluster' with a single cutoff that results in the
    longest processing time, plus some small overhead linear in that processing
    time.
}
\value{
    Returns a data frame. Besides a variable giving compound ID, each of the
    other variables in the data frame will either give the cluster IDs of
    compounds under some clustering cutoff, or the size of clusters that the
    compounds belong to. When N cutoffs are given, in total 2*N+1 variables
    will be generated, with N of them giving the cluster ID of each compound
    under each of the N cutoffs, and the other N of them giving the cluster
    size under each of the N cutoffs. The rows are sorted by the cluster sizes.
}
\author{Y. Eddie Cao, Li-Chang Cheng}
\seealso{\code{\link{cmp.parse1}}, \code{\link{cmp.parse}}, \code{\link{cmp.search}}, \code{\link{cmp.similarity}}}
\examples{
# load sample database from web
db <- cmp.parse("http://bioweb.ucr.edu/ChemMineV2/static/example_db.sdf")
# cluster it
clusters <- cmp.cluster(db, cutoff=0.65)
# cluster using multiple cutoffs
clusters <- cmp.cluster(db, cutoff=c(0.5, 0.85))
# or save the distance before clustering:
clusters <- cmp.cluster(db, cutoff=0.65, save.distances="distmat.rda")
# later, you can load the matrix and pass it to do clustering. Load will load
# the variable 'distmat' that contains the distance matrix
load("distmat.rda")
clusters <- cmp.cluster(db, cutoff=0.60, use.distances=distmat)
}
\keyword{utilities}
