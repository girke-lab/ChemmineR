\name{initDb}
\alias{initDb}
\title{
	Iinitialize SQL Database
}
\description{
	This will ensure that the database connection given is ready for use. If it does not find the tables it needs, it
	will try to create them.
}
\usage{
initDb(handle)
}
\arguments{
  \item{handle}{
	  This can be either a filename, in which case we assume it is the name of an SQLite database and use RSQLite to
	  connect to it, or else any DBI Connection.
}
}
\value{
   Returns a connection object that can be used with other database oriented functions.
}
\author{
   Kevin Horan
}

\seealso{
   RSQLite
}
\examples{
   #create and initialize a new SQLite database
   conn = initDb("test.db")
}
