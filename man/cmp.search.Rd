\name{cmp.search}
\alias{cmp.search}
\title{Search a descriptor database for compounds similar to query compound}
\description{
	Given descriptor of a query compound and a database of compound
	descriptors, search for compounds that are similar to the query
	compound. User can limit the output by supplying a cutoff similarity
	score or a cutoff that limits the number of returned compounds. The
	function can also return the scores together with the compounds.
}
\usage{
    cmp.search(db, query, type=1, cutoff = 0.5, return.score = FALSE, quiet = FALSE,
		    mode = 1, visualize=FALSE, visualize.browse=TRUE, visualize.query=NULL)
}
\arguments{
	\item{db}{The compound descriptor database returned by 'cmp.parse'.}
	\item{query}{The query descriptor, which is usually returned by 
        'cmp.parse1'.}
	\item{type}{Returns results in form of position indices (type=1), 
        named vector with compound IDs (type=2) or data frame (type=3).}
	\item{cutoff}{The cutoff similarity (when cutoff <= 1) or the number of
	maximum compounds to be returned (when cutoff > 1).}
	\item{return.score}{Whether to return similarity scores. If set to
	TRUE, a data frame will be returned; otherwise, only the compounds'
	indices in the database will be returned in the order of decreasing
	scores.}
	\item{quiet}{Whether to disable progress information.}
	\item{mode}{Mode used when computing similarity scores. This value is
	passed to 'cmp.similarity'.}
	\item{visualize}{Whether to visualize the search result in a webpage.}
	\item{visualize.browse}{Whether to open the browser automatically if you choose to visualize the search result.}
	\item{visualize.query}{Filename/URL or a character string containing the  SDF of the query structure if you also want to visualize the query in the search result visualization webpage.}
}
\details{
	'cmp.search' will go through all the compound descriptors in the
	database and calculate the similarity between the query compound and
	compounds in the database. When cutoff similarity score is set,
	compounds having a similarity score higher than the cutoff will be
	returned. When maximum number of compounds to return is set to N via
	'cutoff', the compounds having the highest N similarity scores will be
	returned.

	If 'visualize' is set to a TRUE value, \code{\link{sdf.visualize}} will
	be called to send the search results and the scores to ChemMine
	website. If 'visualize.browse' is set to a TRUE value, the browser will
	open to show the structures in the search result with their
	corresponding scores. Otherwise, a URL pointing to that webpage will be
	printed. By default, 'visualize.query' is not set, and the query
	structure will not be uploaded. If you want that to be included in the
	visualization webpage as well, you must set this argument to a
	character string containing the SDF of the query, or a filename
	pointing to a file containing the SDF of the query. If the character
	string or the file containing multiple SDFs, only the first will be
	considered as the SDF of the query.
}
\value{
  When 'return.score' is set to FALSE, a vector of matching compounds' indices
  in the database will be returned. Otherwise, a data frame will be returned:
  \item{ids}{The indices of matching compounds in the database.}
  \item{scores}{The similarity scores between the matching compounds and the
  query compound}
}
\references{Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", in \emph{J Chem Inf Comput Sci}.}
\author{Y. Eddie Cao, Li-Chang Cheng}
\seealso{\code{\link{cmp.parse1}}, \code{\link{cmp.parse}},
	\code{\link{cmp.search}}, \code{\link{cmp.cluster}},
	\code{\link{cmp.similarity}}, \code{\link{sdf.visualize}}}
\examples{
## Load sample SD file
# data(sdfsample); sdfset <- sdfsample

## Generate atom pair descriptor database for searching
# apset <- sdf2ap(sdfset) 

## Loads same atom pair sample data set provided by library
data(apset) 
db <- apset
query <- db[1]

## Ooptinally, save the db for future use
save(db, file="db.rda", compress=TRUE)

## Search for similar compounds using similarity cutoff
cmp.search(db, query, cutoff=0.2, type=1) # returns index
cmp.search(db, query, cutoff=0.2, type=2) # returns named vector
cmp.search(db, query, cutoff=0.2, type=3) # returns data frame

# you may visualize the search result in ChemMine
\dontrun{cmp.search(db, query, cutoff=10, visualize=TRUE, visualize.browse=FALSE, visualize.query=url)}

## in the next session, you may use load a saved db and do the search:
load("db.rda")
cmp.search(db, query, cutoff=3)
## you may also use the loaded db to do clustering:
cmp.cluster(db, cutoff=0.35)
}
\keyword{utilities}
