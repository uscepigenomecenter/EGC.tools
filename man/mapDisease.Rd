\name{mapDisease}
\alias{mapDisease}
\alias{mapTissue}
\alias{mapHistology}
\alias{mapBatch}

\title{
Map the given information to the information present in DANEUS
}
\description{
Map the given informaion to the information present in the corresponding table in the Epigenome Center Methylation database DANEUS. See Details
}
\usage{
mapDisease(con = NULL, disease)
mapTissue(con = NULL, tissue)
mapHistology(con = NULL, mappings)
mapBatch(con = NULL, mappings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{disease}{
	A character vector of valid TCGA disease types
}
  \item{tissue}{
	A character vector of valid TCGA tissue types
}
  \item{mappings}{
	A data frame containing mapping information like sample names, batch, tissue, disease, histology, etc
}
}
\details{
The functions return a mapping of the provided information to their assigned IDs in DANEUS and other specific information for each function. If the given information is not found in the database the functions try to insert it into the database. If the insertion is unsuccessful the functions return the error message returned by the MySQL INSERT query.

mapDisease and mapTissue take as input a list of disease types and tissue types respectively and map them to the disease and tissue types present in the DISEASE and TISSUE tables in DANEUS and return the corresponding IDs if found or they try to insert the new disease and tissue types into DANEUS first and then return the mapped IDs.

mapHistology takes as input a data frame containing the mapping information and from it first extracts and maps the disease and tissue information and then combines it with the given histologies and maps them to the HISTOLOGY table to retrieve the mapped IDs. If any of the histologies are not present in the table the function tries to add them into DANEUS first and then returns the mapped IDs.

mapBatch takes as input a data frame containing the mapping information and from it first extracts and maps the disease information and then combines it with the given batches and maps them to the BATCH table to retrieve the mapped IDs. If any of the batches are not present in the table the function tries to add them into DANEUS first and then returns the mapped IDs. The function also automatically determines and sets the ordering of any new batches based on the ordering information present in the BATCH table.
}
\value{
An integer vector of mapped IDs
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
\code{\link{dbConnect}}, \code{\link{getDisease}}, \code{\link{insertDisease}}, \code{\link{getTissue}}, \code{\link{insertTissue}}, \code{\link{getHistology}}, \code{\link{insertHistology}}, \code{\link{getBatch}}, \code{\link{insertBatch}}
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
