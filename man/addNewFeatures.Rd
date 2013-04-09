\name{addNewFeatures}
\alias{addNewFeatures}
\title{
   Add New Features
}
\description{
   Adds new features to a database without adding any
   data. Note that if you are loading new data anyway, it
   is much more efficient to use the loadSdf function and
   include the new features then. This function will have
   to read all compounds out of the database first.
}
\usage{
addNewFeatures(conn, featureGenerator)
}
\arguments{
  \item{conn}{
	  A database connection object, such as is returned by \code{\link{initDb}}.
   }
  \item{featureGenerator}{
      A function which returns a data frame containing the
      new features. It may also contain features which are
      already in the database, these will simply be
      ignored. See the description of \code{fct} in
      \code{\link{loadSdf}} for details.
   }
}
\value{
   No value is returned.
}
\author{
   Kevin Horan
}

\seealso{
   \code{\link{loadSdf}} 
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test.db")

	data(sdfsample)

	#just load the data with no features or descriptors
	ids=loadSdf(conn,sdfsample)
   addNewFeatures(conn, function(sdfset) 
					data.frame(MW = MW(sdfset),  
               rings(sdfset,type="count",upper=6, arom=TRUE))
			 )

	unlink("test.db")

}
