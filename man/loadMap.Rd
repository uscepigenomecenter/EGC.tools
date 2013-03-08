\name{loadMap}
\alias{loadMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a mapping file for a tumor from DANEUS
}
\description{
Retrieves information for all the samples for a particular tumor type from DANEUS and generates a mapping file which is required to build and package a TCGA archive
}
\usage{
loadMap(con, disease, platform = "HumanMethylation450")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
	A valid connection object as produced by dbConnect() which provides a connection to DANEUS
}
  \item{disease}{
	Character which specifies the TCGA tumor type for which a mapping file is to be generated
}
  \item{platform}{
	Character specifying the platform for which the data is to be retrieved
}
}
\details{
The \code{loadMap} function sends a \emph{SELECT} query to retrieve all the samples and the other necessary information from DANEUS and generates a valid mapping file which is required as input to \code{packageArchive} to generate a TCGA archive. The \emph{SELECT} query is designed to retrieve data for only those samples that have been completely processed using the specified \code{platform}.
}
\value{
A data frame containing the following information
\item{barcode }{Character vector containing the Illumina barcode for the samples}
\item{uuid }{Character vector containing Universally Unique Identifiers \acronym{UUID} assigned to each sample and were retrieved from the \bold{Biospecimen Metadata Browser}}
\item{TCGA.ID }{Character vector containing the TCGA identifier assigned to each sample}
\item{TCGA.BATCH }{Character vector containing the TCGA Batch associated with each sample}
\item{BATCH.ID }{Interger vector containing the Batch order associated with each TCGA Batch on the TCGA Data Portal}
\item{histology }{Character vector containing the Histology assigned to each sample by the TCGA}
\item{tissue }{Character vector containing the Tissue type from which the sample was derived}
\item{diseaseabr }{Character vector containing the standard Disease Abbreviation assigned to the Tumor type by the TCGA}
\item{plate }{Integer vector containing the Plate Number on which the sample was assayed. This information is useful for assessing batch effects}
\item{well }{Character vector containing the Well position on the Plate where the sample was assayed}
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
\code{\link{dbConnect}}
}
\examples{
\dontrun{
	#Create a connection to DANEUS
	con <-  dbConnect(MySQL(), host, user, pass, db='daneus')
	#Generate a mapping file for the tumor type STAD for the 450k platform
	map <- loadMap(con, "STAD")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
