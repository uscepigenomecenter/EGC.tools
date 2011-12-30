getMetadata <- function(x) 
{ # {{{ skip AtomicLists in elementMetadata(x) DataFrame
  message('Warning: this probably will not work out exactly as you expect...')
  emd = names(elementMetadata(x))
  skip = sapply(emd, function(y) is(elementMetadata(x)[[y]], 'AtomicList'))
  as.data.frame(lapply(emd[-skip], function(z) unlist(elementMetadata(x)[[z]])))
} # }}}

GR2df <- function(GR, keepColumns=FALSE, ignoreStrand=FALSE) 
{ # {{{
  stopifnot(class(GR) == 'GRanges')
  tmp = data.frame(chr=as.vector(seqnames(GR)), start=start(GR), end=end(GR))
  if(!ignoreStrand) tmp = cbind(tmp, strand=as.vector(strand(GR)))
  if(keepColumns) tmp = cbind(tmp, getMetadata(GR))
  return(tmp)
} # }}}
