\name{justSign}
\alias{justSign}
\title{
  Generate MD5sums for a TCGA tumor archive.
}
\description{
  This is a (much) faster version of packageAndSign which is sufficient for 
  running the validator with the -bypass -noremote arguments (i.e. no tarballs).
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
