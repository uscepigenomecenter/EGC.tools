getBarcodes <- function(UUIDs) {
  wsurl <- 'https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/'
  .getBarcode <- function(uuid) fromJSON(getURL(paste0(wsurl, uuid)))['barcode']
  unname(sapply(UUIDs, .getBarcode))
}
