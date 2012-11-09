\name{getDisease}
\alias{getDisease}
\alias{getTissue}
\alias{getHistology}
\alias{getBatch}
\alias{getStatus}
\alias{getPlatform}
\alias{getProject}

\title{
Retrieve information from the various Tables in DANEUS
}
\description{
Retrieve information from the various Tables in the Epigenome Center Methyation database DANEUS. See Details
}
\usage{
getDisease(con = NULL)
getTissue(con=NULL)
getHistology(con=NULL, tissue, disease)
getBatch(con=NULL, disease)
getStatus(con=NULL, status)
getPlatform(con=NULL, platform)
getProject(con=NULL, project)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{tissue}{
	Integer vector of Tissue Type IDs
}
  \item{disease}{
	Integer vector of Disease Type IDs
}
  \item{status}{
	Character vector specifying a Sample's status. All the values should be identical.
}
  \item{platform}{
	Character vector specifying a supported Platform type. All the values should be identical.
}
  \item{project}{
	Character vector specifying the type of Project. All the values should be identical.
}
}
\details{
All the functions listed above send a SELECT query to retrieve relevant information from a particular table in DANEUS.

getDisease and getTissue retrieve all the entries from the DISEASE and TISSUE tables respectively. Currently all the disease types supported are valid acronyms of the various tumor types being investigated by The Cancer Genome Atlas \acronym{TCGA}. This information can be found at \url{https://tcga-data.nci.nih.gov/tcga}.

getHistology retrieves all the entries from the HISTOLOGY table for only a given set of disease and tissue types.

getBatch retrieves all the entries from the BATCH table for only a given set of disease types.

getStatus, getPlatform and getProject retrieve the corresponding ID for the queried status, platform and project from the STATUS, PLATFORM and PROJECT tables. Some more information about these 3 particular tables is provided below.

The STATUS table contains only four valid statuses as of now. They are \sQuote{BS Converted, QC Failed}, \sQuote{Nanodrop Pass}, \sQuote{Methylation QC Pass}, \sQuote{Infinium}.

Currently only two platforms are supported in the PLATFORM table and they are \sQuote{HumanMethylation27}, \sQuote{HumanMethylation450}.

The Methylation database currently only contains sample from the following two PROJECT types which are \sQuote{TCGA} and \sQuote{EGC}.
}
\value{
getDisease and getTissue return a data fram containing the following information:
\item{name }{Character vector containing the disease or tissue types}
\item{id }{Integer vector containing the IDs assigned to the correponding disease and tissue types in \option{name}}



getHistology returns the following data frame:
\item{name }{Character vector of Histology types}
\item{id }{Integer vector of IDs assigned to the histology types}
\item{tissue }{Integer vector of Tissue IDs}
\item{disease }{Integer vector of Disease IDs}



getBacth returns the following data frame:
\item{disease }{Integer vector of Disease IDs}
\item{id }{Integer vector of IDs assigned to the batches}
\item{ordering }{Integer vector containing the batch number assigned to a particular batch on the TCGA Data Portal}
\item{batch }{Character vector containing the actual TCGA Batch assigned to a Batch in the BCR shipping manifest}



getStatus, getPlatform and getProject each return an Integer vector containing the ID assigned to the queried status, platform or project respectively
}
\references{
The Cancer Genome Atlas Project \url{https://tcga-data.nci.nih.gov/tcga/}
}
\author{
Moiz Bootwalla
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbConnect}}, \code{\link{dbGetQuery}}
}
\examples{
\dontrun{
	#Create a connection to DANEUS
	con <-  dbConnect(MySQL(), host, user, pass, db='daneus')
	disease <- getDisease(con)
	tissue <- getTissue(con)
	histology <- getHistology(con, disease=disease, tissue=tissue)
	batch <- getBatch(con, disease=disease)
	getStatus(con, status="Infinium")
	getPlatform(con, platform="HumanMethylation450")
	getProject(con, project="TCGA")
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
