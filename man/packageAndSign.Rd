\name{packageAndSign}
\title{
  Generate MD5sums and tarballs for a TCGA tumor archive.
}
\description{
  Read the source code, or at least the arguments.  This is the last function
  typically called before uploading an archive via SFTP. 
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
