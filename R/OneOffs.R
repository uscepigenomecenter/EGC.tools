## Various functions to build TCGA packages and such
loadOneOff <- function(mapping, label, platform='HumanMethylation450', path='/auto/uec-02/shared/production/methylation', project='tcga') {
  # {{{
  require(methylumi)
  platform.path = paste(path, 
                        ifelse(grepl('HumanMethylation450',platform,ignore=T),
                               'meth450k', 'meth27k'), 'raw', sep='/')
  oneoff.dir = paste(gsub('raw', 'processed', platform.path), project, label, sep='/')
  save.dir = gsub('raw', 'MethyLumiSets', platform.path) 
  if(!file.exists(oneoff.dir)) dir.create(oneoff.dir)
  stopifnot('barcode' %in% names(mapping))
  rownames(mapping) = mapping$barcode
  for(b in as.character(mapping$barcode)) {
    bdir = paste(platform.path, substr(b, 1, 10), b, sep='/')
    glob = paste(bdir, '*.idat', sep='')
    system(paste('ln', glob, oneoff.dir))
    IDATs = unlist(list.files(path=oneoff.dir, patt=b))
    print(paste('Linked', paste(IDATs, collapse=' and '), 'to', oneoff.dir))
  }
  setwd(oneoff.dir)
  methylumIDAT(mapping$barcode, n.sd=T, oob=T, parallel=T)
} # }}}

writeOneOffExcel <- function(x, label) {
  # {{{
  message("Dumping betas, M, U, p-values and controls to Excel sheets...")
  stop("This function ain't implemented yet")
} # }}}

writeOneOffCSVs <- function(x, label) {
  # {{{
  message("Dumping betas, M, U, and p-values, and controls to CSV files...")
  stop("This function ain't implemented yet")
} # }}}

plotOneOffQC <- function(x, label) {
  # {{{
  message("Generating QC plots...")
  stop("This function ain't implemented yet")
} # }}}

saveOneOff<- function(x, label) {
  # {{{
  stop("This function ain't implemented yet")
} # }}}

