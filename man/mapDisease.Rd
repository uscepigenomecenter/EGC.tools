\name{mapDisease}
\alias{mapDisease}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve Disease ids from DANEUS
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mapDisease(con = NULL, disease)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
%%     ~~Describe \code{con} here~~
}
  \item{disease}{
%%     ~~Describe \code{disease} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (con = NULL, disease) 
{
    require("RMySQL")
    if (is.null(con)) 
        stop("Please provide a Database connection object: see ?dbConnect")
    disease.db <- getDisease(con)
    map <- match(disease, disease.db[, "name"])
    if (any(is.na(map))) {
        new.disease <- unique(disease[which(is.na(map))])
        new.disease <- ifelse(length(new.disease) > 1, paste(new.disease, 
            collapse = "','"), new.disease)
        message(paste("Inserting '", new.disease, "' into the DISEASE Table", 
            sep = ""))
        if (insertDisease(con, new.disease)) {
            disease.db <- getDisease(con)
            map <- match(disease, disease.db[, "name"])
        }
        else {
            stop(paste("There was an error inserting '", new.disease, 
                "' into the database", sep = " "))
        }
    }
    disease.id <- disease.db[, "id"][map]
    return(disease.id)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
