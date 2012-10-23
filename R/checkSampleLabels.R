checkSampleLabels <- function(map)
{ # {{{
  stopifnot('histology' %in% names(map))
  stopifnot('TCGA.ID' %in% names(map))
  labeling = c('01'='Tumor', '11'='Normal', '20'='Control')
  map$designation = as.factor(labeling[ substr(map$TCGA.ID, 14, 15) ])
  if(length(levels(map$designation)) > 1){
    map$designation = relevel(map$designation, 
                              which(levels(map$designation)=='Tumor'))
  }
  names(map)[which(names(map) == 'designation')] = 'designation.by.TCGA.ID'
  cat("\n","Please check the following table for any inconsistencies:","\n\n")
  with(map, print(table(toupper(histology), designation.by.TCGA.ID)))
} # }}}
