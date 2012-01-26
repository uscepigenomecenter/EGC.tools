\name{mapBatches}
\alias{mapBatches}
\title{
  Given a pile of Dan's Excel files, attempt to merge them into a map.
}
\description{
  If the function squawks, do what it says, and fix the offending elements.
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
