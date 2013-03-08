\name{linkRawData}
\alias{linkRawData}
\title{
  From a disease map (a data.frame), link up all of the raw IDAT and SDF files.
}
\description{
  By default, this creates $HOME/meth450k/raw/DISEASEABR/ and fills it up. The function first checks if all the barcodes are present in $HOME/meth450k/raw and throws an error listing the missing barcodes if there are any
}
\usage{
linkRawData(map, unlink.old.files = TRUE, platform = "HumanMethylation450", path = "/export/uec-gs1/laird/shared/production/methylation", logfile = NULL, idatPath = NULL)
}
\arguments{
  \item{map}{Data frame containing at the minimum the element \code{barcode} which contains the Illumina barcodes of the samples}
  \item{unlink.old.files}{Logical specifying whether or not existing directories should be erased or not}
  \item{platform}{Character specifying the platform for which the IDAT files have to be linked}
  \item{path}{Character specifying the default directory where to look for the raw data}
  \item{logfile}{Character specifying the name of the log file to which to write out the process log. If it is \code{NULL} then no log file is generated and the process log is sent to \env{STDOUT}}
  \item{idatPath}{Character specifying the directory to look for the IDATs in if they are not in the standard directory structure. See \emph{Warning} for details}
}
\section{Warning }{
  This function is tailor made for the directory structure implemented on USC's High Performance Compute Cluster \acronym{HPCC}. See \code{\link{buildArchive}} function for more details regarding the organisation of the directory structure
}
\author{
Tim Triche Jr. and Moiz Bootwalla
}

