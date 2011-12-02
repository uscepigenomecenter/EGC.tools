# Write processing logs.

getTidyHistory <- function(x) 
{ # {{{ 
  h = getHistory(x)
  th = h[,c('submitted','command')]
  th$submitted = sapply(as.character(th$submitted), function(x) {
                          strsplit(x,' ')[[1]][2]
                        })
  th$command = as.character(th$command)
  return(th)
} # }}}

writeTidyHistory <- function(x, filepath='.') 
{ # {{{ usually will want this in AUX
  disease = unique(x$diseaseabr)[1]
  filename = paste(disease,'processing','log','txt',sep='.')
  write.table(getTidyHistory(x), 
              file=paste(filepath, filename, sep='/'),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
} # }}}

