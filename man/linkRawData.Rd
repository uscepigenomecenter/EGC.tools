\name{linkRawData}
\title{
  From a disease map (a data.frame), link up all of the raw IDAT and SDF files.
}
\description{
  By default, this creates $HOME/meth450k/raw/DISEASEABR/ and fills it up.
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
