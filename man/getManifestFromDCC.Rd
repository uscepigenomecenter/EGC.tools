\name{getManifestFromDCC}
\alias{getManifestFromDCC}

\title{
Fetch MANIFEST files from the DCC
}
\description{
The function is a utility function to retrieve the MANIFEST files for the specified batch and version for a TCGA Tumor type from the Data Portal. This function is helpful while creating revision archives (the revision archive creation functionality is not yet completely implemented)
}
\usage{
getManifestFromDCC(version = "0", level = "1", disease, batch, platform = "HumanMethylation450")
}

\arguments{
  \item{version}{
	Character specifying the version of the archive whose MANIFEST files need to be fetched
}
  \item{level}{
	Character vector specifying the archive levels whose MANIFEST files need to be fetched
}
  \item{disease}{
	Character specifying the TCGA tumor type
}
  \item{batch}{
	Integer vector specifiying the batch numbers for which the MANIFEST files are to be retrieved. Note that these should be the same as that on the TCGA Data Portal and not be the TCGA specified Batch number. This information is generally present in the TUMOR.mappings.csv file in the corresponding aux directory for the tumor on the TCGA Data Portal
}
  \item{platform}{
	The assay platform for which the MANIFEST files are being retrieved
}
}

\author{
Moiz Bootwalla
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
