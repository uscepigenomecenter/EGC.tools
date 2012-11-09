\name{insertDisease}
\alias{insertDisease}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Insert new disease types into DANEUS's DISEASE Table
}
\description{
Insert new disease types and assign them an ID in the Epigenome Center Methylation Database DANEUS's DISEASE Table
}
\usage{
insertDisease(con = NULL, disease)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{disease}{
	A single character value containing the disease types to be inserted separated by the separator "','". For example DISEASE1','DISEASE2','DISEASE3
}
}
\details{
The function tries to insert the given diseases types into the DISEASE Table in DANEUS and if it is unsuccessful it returns the error returned by the MySQL INSERT query
}
\value{
Logical indicating success or failure of MySQL INSERT operation
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

\section{Warning }{
This is a lower level function which requires the list of new diseases to be inserted to be provided in a specific format. Other higher level functions like \code{\link{mapDisease}} which call it do the requisite formatting to the input automatically. Use this function with caution especially when trying to insert multiple new disease types into the DISEASE Table.
}

\seealso{
\code{\link{dbConnect}}, \code{\link{dbSendQuery}}, \code{\link{mapDisease}}
}
\examples{
\dontrun{
	#Create a connection to DANEUS
	con <-  dbConnect(MySQL(), host, user, pass, db='daneus')
	disease <- paste("D1", "D2", "D3", sep="','")
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
