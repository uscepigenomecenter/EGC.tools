\name{mapDisease}
\alias{mapDisease}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map a set of Diseases to their assigned IDs in DANEUS
}
\description{
Given a vector of disease types, retrieve and map their assigned IDs from DANEUS
}
\usage{
mapDisease(con = NULL, disease)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{disease}{
	A character vector of valid TCGA disease types
}
}
\details{
The function returns a mapping of the provided diseases to their assigned IDs in DANEUS. If a disease type is not found in the database the function tries to insert it into the database and assign it an ID. If the insertion is unsuccessful the function returns the error message returned by the MySQL INSERT query.
}
\value{
An integer vector of disease IDs
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Moiz Bootwalla
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbConnect}},
}
\examples{
\dontrun{
	disease <- rep(c("LUSC", "THCA", "BLCA"), each=3)
	#Create a connection to DANEUS
	con <-  dbConnect(MySQL(), host, user, pass, db='daneus')
	id <- mapDisease(con, disease)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
