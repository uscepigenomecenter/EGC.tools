# tidy up Dan's mapping spreadsheets (i.e. BeadChip Decoder sheets)

canonicalizeExcelFile <- function(filename) {
  # {{{ 
  require('xlsx')
  message(paste("Reading first worksheet in", filename, "as mappings..."))
  map = read.xlsx(filename, 1, stringsAsFactors=FALSE) 
  newmap = canonicalizeMapping(map)
  detach('package:xlsx', unload=TRUE) # Java fucks things up otherwise
  return(newmap)
} # }}}

canonicalizeMapping <- function (map, synonyms = NULL) {
  # {{{
    map = map[which(!is.na(map[,1])),] # blunt but effective
    if (is.null(synonyms)) {
      synonyms = list(barcode = c("Complete.Barcode", "Complete.Barcode.ID"),
                      histology = c("Histology", "X.Histology"), 
                      tissue = c("Tissue.Type"), 
                      diseaseabr = c("Disease.Abbreviation"), 
                      TCGA.BATCH = c("Batch", "TCGAbatch"), 
                      TCGA.ID = c("Biospecimen.Barcode.Side"), 
                      source = c("Plate.Well"))
    }
    for (i in setdiff(names(synonyms), names(map))) {
      for (j in synonyms[[i]]) if (j %in% names(map)) {
        map[[i]] = map[[j]]
      }
    }
    separators = list(barcode = "_", source = " ")
    components = list(barcode = c("Sentrix_ID", "Terminus"),
                      source = c(PLATE = "Plate", WELL = "Well.Position"))
    for (k in setdiff(names(components), names(map))) {
      if (all(components[[k]] %in% names(map))) {
        piece.names = names(components[[k]])
        if (!is.null(piece.names)) {
          map[[k]] = paste(c(rep(piece.names[[1]], dim(map)[1])), 
                           map[[components[[k]][1]]], 
                           c(rep(piece.names[[2]], dim(map)[1])), 
                           map[[components[[k]][2]]], 
                           sep = separators[[k]])
        } else {
          map[[k]] = paste(map[[components[[k]][1]]], map[[components[[k]][2]]],
                           sep = separators[[k]])
        }
      } else {
        message("The minimum subset of columns required is:")
        message("Biospecimen.Barcode.Side, Histology, Disease.Abbreviation,")
        message("Sentrix_ID, Terminus (to successfully map tumors and normals)")
        message(paste("Could not assemble column",k,"from component pieces..."))
      }
    }
    map$name = map$TCGA.ID # ugly horrible disgusting hack
    if( !all( names(synonyms) %in% names(map) ) ) {
      message(paste('Missing',
                    paste( setdiff(names(synonyms), names(map)), collapse=', '),
                    'from map'))
      message('Possible synonyms:')
      for( nm in setdiff(names(synonyms), names(map)) ) {
        message(paste(nm, paste(synonyms[[nm]], collapse=' or ')))
      }
    } else { 
      return(map[which(!is.na(map$barcode)), names(synonyms)])
    }
} # }}}
