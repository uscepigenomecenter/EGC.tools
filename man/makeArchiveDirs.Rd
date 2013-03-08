\name{makeArchiveDirs}
\alias{makeArchiveDirs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low-level function used to create aux, Levels 1, 2, 3 and magetab directories for the TCGA archive being packaged
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeArchiveDirs(map, version = "0", base = NULL, magetab.version = NULL, lvls = c("aux", "Level_1", "Level_2", "Level_3", "mage-tab"), platform = "HumanMethylation450")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
	A MethyLumiSet object or a data frame containing metadata about the TCGA archive being built
}
  \item{version}{
	Character specifying the version of the archive being generated
}
  \item{base}{
	The path to the base directory where the archive will be created
}
  \item{magetab.version}{
	The version of the mage-tab directory to create. If \code{NULL} it defaults to \code{version}
}
  \item{lvls}{
	Character vector containing the levels for which to create directories
}
  \item{platform}{
	The assay platform for which the archive is being created
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List containing the paths to the aux, Level 1,2,3 and mage-tab directories
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tim Triche Jr. and Moiz Bootwalla
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
