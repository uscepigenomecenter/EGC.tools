require(EGC.tools) 

test_getBarcodesAndUUIDs <- function() {
  barcode <- 'TCGA-AB-2977'
  UUID <- 'f95063bc-a640-4980-99a7-84645902daf4'
  checkTrue( getBarcodes(getUUID(barcode)) == getBarcodes(UUID) )
  checkTrue( getUUIDs(getBarcode(UUID)) == getUUIDs(barcode) )
}
