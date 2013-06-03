\name{insertDisease}
\alias{insertDisease}
\alias{insertTissue}
\alias{insertHistology}
\alias{insertBatch}
\alias{insertSample}

\title{
Insert information into the various tables in DANEUS
}
\description{
Insert new information into various tables in the Epigenome Center Methylation Database DANEUS. See Details
}
\usage{
insertDisease(con = NULL, disease)
insertTissue(con = NULL, tissue)
insertHistology(con = NULL, histology)
insertBatch(con = NULL, batch)
insertSample(con = NULL, mappings, platform="HumanMethylation450", project="TCGA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{disease}{
	A single character value containing the disease types to be inserted separated by the separator "','". For example DISEASE1','DISEASE2','DISEASE3
}
  \item{tissue}{
	A single character value containing the tissue types to be inserted separated by the separator "','". For example TISSUE1','TISSUE2','TISSUE3
}
  \item{histology}{
	A character vector where each element is of the form "'HISTOLOGY','TISSUE_ID','DISEASE_ID'". The TISSUE_ID and DISEAE_ID are retrieved from DANEUS
}
  \item{batch}{
	A character vector where each element is of the form "'DISEASE_ID','BATCH','ORDERING'". The DISEAE_ID is retrieved from DANEUS
}
  \item{mappings}{
	A data frame containing columns that are required to populate the SAMPLE table in DANEUS. See details
}
  \item{platform}{
	Character specifying a supported Platform type. Set to 'HumanMethylation450' by default
}
  \item{project}{
	Character specifying the type of Project. Set to 'TCGA' by default
}
}
\details{
All the functions listed above send an INSERT query to insert the specified information into a particular table in DANEUS. The functions try to insert the specified information into DANEUS and if they are unsuccessful they return the error returned by the MySQL INSERT query.

insertDisease and insertTissue insert the provided disease and tissue types into the DISEASE and TISSUE tables respectively.

insertHistology inserts the histology for a given set of disease and tissue types. The disease and tissue ids from the respective DISEASE and TISSUE tables and not their actual names are required for the INSERT query to work successfully.

insertBatch inserts the specified batches into the BATCH table for a particular disease type. The disease id from the DISEASE table is required for the INSERT query to work successfully.

insertSample inserts the specified samples into the SAMPLE table. The function requires a data frame containing minimal sample information like the TCGA.ID and shipping date to work successfully. Other columns in the SAMPLE table like barcode, uuid, plate and well are optional and if not provided are set to NULL and must be updated later on to generate a functional mapping file using \code{\link{loadMap}}. The status if not specified is set to 'Nanodrop Pass'. Refer to the STATUS table to check out currently accepted statuses.
}
\value{
Logical indicating success or failure of MySQL INSERT operation. If query fails then the error produced by the MySQL INSERT statement is also returned.
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
These are lower level functions which require the input values to be provided in a specific format. Other higher level functions like \code{\link{mapDisease}} call their respective insert functions like \code{\link{insertDisease}} and do the requisite formatting to the input automatically. Use these functions with caution especially when trying to insert multiple values into any of the tables.
}

\seealso{
\code{\link{dbConnect}}, \code{\link{dbSendQuery}}, \code{\link{dbWriteTable}}, \code{\link{mapDisease}}, \code{\link{mapTissue}}, \code{\link{mapHistology}}, \code{\link{mapBatch}}
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
