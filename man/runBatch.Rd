\name{runBatchByID}
\alias{runBatchByID}
\alias{runBatchByName}

\title{
Legacy function to process an archive one batch at a time
}
\description{
This function was initially written to process BRCA archives. Improved functionality in \code{buildArchive} should now be used instead of these legacy functions. Refer to the User Manual for details on processing large datasets.
}
\usage{
runBatchByID(map, batch.id, base = NULL, platform = "HumanMethylation450")
runBatchByName(map, name, base = NULL, platform = "HumanMethylation450")
}

\author{
Tim Triche Jr.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
