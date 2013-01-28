loadOneOff <- function(mapping, label, platform='HumanMethylation450', path='/export/uec-gs1/laird/shared/production/methylation', project='tcga') 
{ # {{{
  require(methylumi)
  platform.path = paste(path, 
                        ifelse(grepl('HumanMethylation450',platform,ignore=T),
                               'meth450k', 'meth27k'), 'raw', sep='/')
  oneoff.dir = paste(gsub('raw', 'processed', platform.path), project, label, sep='/')
  save.dir = gsub('raw', 'MethyLumiSets', platform.path) 
  if(!file.exists(oneoff.dir)) dir.create(oneoff.dir)
  stopifnot('barcode' %in% names(mapping))
  rownames(mapping) = mapping$barcode
  checkDirs <- function(b){
	  return(file.exists(b))
  }

  check <- sapply(substr(mapping$barcode,1,10), function(x){checkDirs(paste(platform.path, x, sep='/'))})
  if(!all(check)){
	  failed <- paste(unique(names(check[which(check == FALSE)])), collapse=",")
	  print(paste("Barcodes", failed, "were not found"))
	  return(NULL)
  }else{
	  for(b in as.character(mapping$barcode)) {
		  bdir = paste(platform.path, substr(b, 1, 10), b, sep='/')
		  glob = paste(bdir, '*.idat', sep='')
		  system(paste('ln', glob, oneoff.dir))
		  IDATs = unlist(list.files(path=oneoff.dir, patt=b))
		  print(paste('Linked', paste(IDATs, collapse=' and '), 'to', oneoff.dir))
	  }
	  setwd(oneoff.dir)
	  methylumIDAT(mapping$barcode, n.sd=T, oob=T, parallel=T)
  }
} # }}}

runOneOff <- function(label, rootdir=NULL, mapdir=NULL, jobdir=NULL) 
{ # {{{

  ## {{{ default rootdir, mapdir, and jobdir settings
  if(is.null(rootdir)) 
    rootdir <- "/export/uec-gs1/laird/shared/production/methylation/meth450k"
  if(is.null(mapdir)) 
    mapdir <- paste0(rootdir, '/mappings') 
  if(is.null(jobdir)) 
    jobdir <- paste0(rootdir, '/processed/internal/', label)
  ## }}}

  ## read in the chip mappings and any additional pData
  mapfile <- paste0(mapdir, '/', label, '.mapping.csv')
  mapping <- read.csv(mapfile)
  stopifnot(nrow(mapping) > 0 && ncol(mapping) > 0)

  ## load the data from IDAT files
  methyldata <- loadOneOff(mapping, label, project='internal')
  stopifnot(class(methyldata) %in% c('MethyLumiSet','MethylSet'))

  ## save the raw dataset
  rawrda <- paste0(jobdir, '/', label, '.raw.rda')
  save(methyldata, file=rawrda)

  ## write out the raw beta values
  rawbetas <- paste0(jobdir, '/', label, '.betas.raw.txt')
  write.table(betas(methyldata), sep="\t", file=rawbetas)

  ## preprocess: background correct and dye-bias equalize
  methyldata <- stripMethyLumiSet(methylumi.bgcorr(methyldata, 'noob'))
  methyldata <- normalizeMethyLumiSet(methyldata) ## should we set a referent?

  ## save the processed dataset
  fixedrda <- sub('raw', 'corrected', rawrda) ## just swap words :-)
  save(methyldata, file=fixedrda)

  ## write out the corrected beta values
  fixedbetas <- sub('raw', 'corrected', rawbetas) ## as above
  write.table(betas(methyldata), sep="\t", file=fixedbetas)

  ## write out the corrected M intensities
  fixedMintensities <- sub('betas', 'methylated', fixedbetas)
  write.table(methylated(methyldata), sep="\t", file=fixedMintensities)

  ## write out the corrected U intensities
  fixedUintensities <- sub('betas', 'unmethylated', fixedbetas)
  write.table(unmethylated(methyldata), sep="\t", file=fixedUintensities)

  ## write out the p values (which don't change!)
  fixedPvals <- sub('betas', 'pvalues', fixedbetas)
  write.table(pvals(methyldata), sep="\t", file=fixedPvals)

  ## write out the sequence of commands performed for this job
  processinghistory <- paste0(jobdir, '/', label, '.history.R')
  savehistory(processinghistory)

  ## if we got here:
  message(paste('Processing completed for', label))
  message(paste('Results may be found in', jobdir))

} # }}}
