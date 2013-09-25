getUUIDs <- function(barcodes) {
  require(RCurl)
  require(RJSONIO)
  wsurl <- 'https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/barcode/'
  .getUUID <- function(barcode) {
    fromJSON(getURL(paste0(wsurl, barcode)))[['uuidMapping']][['uuid']]
  }
  unname(sapply(barcodes, .getUUID))
}
