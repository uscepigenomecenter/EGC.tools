\name{buildArchive}
\alias{buildArchive}
\alias{EGC.tools}
\title{
  Build a TCGA archive with aux, level 1, level 2, level 3, and mage-tab dirs. 
}
\description{
  Read the source code.  This can be optimized by splitting up the map file and
  then loading in each batch separately, although most of the time I used a moby
  combined dataset to run it.  That is not strictly necessary, however.

  The function will create dummy tarballs and add manifests to each directory
  after writing out the levels, but you should still check on the aux directory
  to ensure that everything has been included, prior to uploading.
}
\examples{

## map <- loadMap(KIRC) # or
## map <- mapBatches(list.files(patt='xls'))
## linkRawData(KIRC)
## setwd('~/meth450k/raw/KIRC') # from above
## KIRC <- methylumIDAT(map, parallel=TRUE)
## table(KIRC$histology)
## 
## # move any old archives out of the way first...
##
## buildArchive(KIRC, version='2')
## packageAndSign(KIRC)

}
