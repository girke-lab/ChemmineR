\name{getCompoundNames}
\alias{getCompoundNames}
\title{
	Get Compound Names
}
\description{
	Fetch the names of the given compound ids, if they exist
}
\usage{
getCompoundNames(conn, compoundIds)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{compoundIds}{
     A vector of compound ids.
   }
}
\value{
   Returns a vector of compound names.The rownames will be the compound ids. Compound ids
   not found, or for which a name is not defined, will be represented as NA.
}
\author{
   Kevin Horan
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test2.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)

   getCompoundNames(conn,ids[1:3])
	unlink("test3.db")

}
