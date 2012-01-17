\name{loadMap}
\alias{loadMap}
\title{
  Load the existing mapping file for a disease and canonicalize it if needed.
}
\description{
  Read the source code, or at least the arguments.
}
\examples{

## map <- loadMap(KIRC) # or
## map <- mapBatches(list.files(patt='xls'))
## linkRawData(KIRC)
## setwd('~/meth450k/raw/KIRC')
## KIRC <- methylumIDAT(map, parallel=TRUE)
## table(KIRC$histology)
## buildArchive(KIRC)
## packageAndSign(KIRC)

}
