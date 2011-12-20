## FIXME: the new world order
##
## Level 0 / AUX: processing logs, scripts, .SDF files, mappings
## Level 1: IDAT files ( two lines per sample in the SDRF )
## Level 2: background-corrected probe intensities ( M, U, Pval )
## Level 3: SNP10-masked beta values ( ADF points to GEO platform descriptor )
##
reloadEgcTools <- function()
{ # {{{
  base::detach(package:EGC.tools, unload=TRUE)
  library(EGC.tools)
} # }}}

## load and canonicalize
loadMap <- function(disease) 
{ # {{{ 
  map = read.csv( paste(disease, 'mappings', 'csv', sep='.'), 
                  stringsAsFactors=F,
                  row.names=1)
  map = canonicalizeMapping(map)
  map$BATCH.ID = as.character(as.numeric(as.factor(map$TCGA.BATCH)))
  return(map)
} # }}}

## process a tumor's worth of Excel files from Dan
mapBatches <- function(xls.files, parallel=FALSE, link.raw=FALSE) 
{ # {{{
  if(parallel) {
    require(parallel)
    map = do.call(rbind, mclapply(xls.files, canonicalizeExcelFile))
  } else { 
    map = do.call(rbind, lapply(xls.files, canonicalizeExcelFile))
  }
  map$BATCH.ID = as.numeric(as.factor(map$TCGA.BATCH))
  if(link.raw == TRUE) linkRawData(map, unlink.old.files=FALSE)
  return(map)
} # }}}

## not currently using this, since run-by-batch seems to be less trouble
loadBatches <- function(map) ## FIXME: use ff for this step, or per-batch sets
{ # {{{
  require(methylumi)
  if(nrow(map) > 97) { # don't want to have single-chip overhangs
    message('Splitting into 96-chip batches... this WILL thrash with low RAM!')
    map.extra = map[ -c(1:96), ] 
    map = map[ c(1:96), ]
    ##
    ## FIXME: use temporary files or SOMETHING to avoid this...
    ##
    return(combine( loadBatches(map), loadBatches(map.extra) ))
  } else { 
    mlumi = methylumIDAT(map, parallel=TRUE)
    mlumi = stripOOB(methylumi.bgcorr(mlumi, 'noob'))
    return( mlumi )
  }
} # }}}

## e.g. runBatchByName(map.UCEC,'137') is the same as runBatchByID(map.UCEC,10)
runBatchByName <- function(map, name, base=NULL, platform='HumanMethylation450')
{ # {{{
  stopifnot('TCGA.ID' %in% names(map))
  batch.id = which(levels(as.factor(map$TCGA.ID)) == as.character(name))
  runBatchByID(map, batch.id, base=base, platform=platform)
} # }}}

## e.g. runBatchByID(map.UCEC, 1) is the same as runBatchByName(map.UCEC, '75')
runBatchByID <- function(map, batch.id,base=NULL,platform='HumanMethylation450')
{ # {{{
  dirs = list()
  oldwd = getwd()
  stopifnot('TCGA.ID' %in% names(map))
  if(!('BATCH.ID' %in% names(map))) { # {{{
    map$BATCH.ID = as.numeric(as.factor(map$TCGA.BATCH))
  } # }}}
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}
  disease = unique(map$diseaseabr)
  stopifnot(length(disease) == 1)
  dirs$disease = paste(base, 'raw', disease, sep='/')
  dirs$archive = paste(base, 'tcga', disease, sep='/')
  diseasestub = paste('jhu-usc.edu', disease, sep='_')
  dirs$processed = paste(base, 'MethyLumiSets', disease, sep='/')
  batchstub = paste('batch', levels(map$TCGA.BATCH)[batch.id], sep='')
  setwd(dirs$processed)
  if(file.exists(paste(batchstub, 'rda', sep='.'))) { # {{{ load it
    load(paste('batch',levels(map$TCGA.BATCH)[batch.id],'.rda',sep=''))
    batch.data = get(batchstub) # cute trick, no? 
    if(!('BATCH.ID' %in% varLabels(batch.data))) {
      batch.data$BATCH.ID = batch.id 
    } # }}}
  } else { # {{{ process and import it
    setwd(dirs$disease)
    batch.map = map[ which(map$BATCH.ID == batch.id), ] 
    batch.data = stripOOB(methylumi.bgcorr(methylumIDAT(batch.map)))
    gc(,T)
  } # }}}
  writeBatch(batch.data, batch.id)
  setwd(oldwd)
} # }}}

## for breaking up big tumors by batch
writeBatch <- function(x,batch.id,version='0',base=NULL,parallel=F,lvls=c(1:3))
{ # {{{ assuming that the full map will be used for aux and mage-tab

  dirs = list()
  disease = unique(x$diseaseabr)
  stopifnot(length(disease) == 1)
  stopifnot('TCGA.BATCH' %in% varLabels(x))
  if(!('BATCH.ID' %in% varLabels(x))) { # {{{
    x$BATCH.ID = as.numeric(as.factor(x$TCGA.BATCH))
  } # }}}
  message(paste('Writing batch ', batch.id, '...',sep=''))
  platform = gsub('k$','',gsub('Illumina','',annotation(x))) 
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}
  diseasestub = paste('jhu-usc.edu', disease, sep='_')
  stopifnot('TCGA.ID' %in% varLabels(x))
  sampleNames(x) <- x$TCGA.ID
  dirs$raw = paste(base, 'raw', disease, sep='/')
  pkged = dirs$archive = paste(base, 'tcga', disease, sep='/')
  message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')

  if(platform == 'HumanMethylation27') { # {{{
    base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # }}}
  } else { # {{{
    base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) 
  } # }}}

  if(1 %in% lvls) { # {{{ level 1: IDAT files
    dirs$level_1 = paste(pkged,paste(diseasestub, platform, 'Level_1', batch.id,
                                     version, '0', sep='.'), sep='/')
    stopifnot('TCGA.BATCH' %in% varLabels(x))
    message(paste('Creating directory', dirs$level_1, '...'))
    if(!file.exists(dirs$level_1)) dir.create(dirs$level_1)
    oldwd = getwd()
    setwd(dirs$level_1)
    batchname = unique(x$TCGA.BATCH[ which(x$BATCH.ID == batch.id) ])
    stopifnot(length(batchname)==1)
    subjects = sampleNames(x)[ which(x$BATCH.ID == batch.id) ]
    for(s in subjects) { # {{{
      xs = which(sampleNames(x) == s)
      message(paste("Copying IDAT files for sample", which(subjects == s),
                    "of", length(subjects), "in TCGA batch", batchname))
      stub = paste(x$barcode[xs], '_', sep='')
      for(i in paste(stub, c('Grn.idat','Red.idat'), sep='')) {
        file.link(paste(dirs$raw, i, sep='/'), paste('.', i, sep='/'))
      }
    } # }}}
    setwd(oldwd)
  } # }}}

  if(2 %in% lvls) { # {{{ level 2: background-corrected M and U intensities
    dirs$level_2 = paste(pkged,paste(diseasestub, platform, 'Level_2', batch.id,
                                     version, '0', sep='.'), sep='/')
    stopifnot('TCGA.BATCH' %in% varLabels(x))
    message(paste('Creating directory', dirs$level_2, '...'))
    if(!file.exists(dirs$level_2)) dir.create(dirs$level_2)
    oldwd = getwd()
    setwd(dirs$level_2)
    batchname = unique(x$TCGA.BATCH[ which(x$BATCH.ID == batch.id) ])
    stopifnot(length(batchname)==1)
    subjects = sampleNames(x)[ which(x$BATCH.ID == batch.id) ]

    write.level2 <- function(s) { # {{{
      xs = which(sampleNames(x) == s)
      message(paste("Writing level 2 data for sample", which(subjects == s),
                    "of", length(subjects), "in TCGA batch", batchname))
      lvl2data = data.frame( M=methylated(x)[,xs], U=unmethylated(x)[,xs] )
      rownames(lvl2data) = featureNames(x)
      dump.file = paste(paste(diseasestub,platform,b,'lvl-2',s,'txt',sep='.'))
      headers1 = paste('Hybridization REF', s, s, sep="\t")
      headers2 = paste('Composite Element REF', 'Methylated_Intensity', 
                       'Unmethylated_Intensity', sep="\t")
      cat(headers1, "\n", sep='', file=dump.file)
      cat(headers2, "\n", sep='', file=dump.file, append=TRUE)
      write.table(lvl2data, file=dump.file, append=TRUE, quote=FALSE,
                            row.names=TRUE, col.names=FALSE, sep="\t")
      return('done')
    } # }}}

    b = batch.id
    if(parallel) {
      if(!require(parallel)) require(multicore)
      results <- mclapply(subjects, write.level2)
    } else { 
      results <- lapply(subjects, write.level2)
    }

    gc(,T)
    setwd(oldwd) 

  } # }}}

  if(3 %in% lvls) { # {{{ level 3: masked beta values and p-values
    dirs$level_3 = paste(pkged,paste(diseasestub, platform, 'Level_3', batch.id,
                                     version, '0', sep='.'), sep='/')
    message(paste('Creating directory', dirs$level_3, '...'))
    if(!file.exists(dirs$level_3)) dir.create(dirs$level_3)
    setwd(dirs$level_3)
    betas(x) <- methylated(x)/total.intensity(x)
    if( !('SNP10' %in% fvarLabels(x)) ) {
      data(SNPs) # can upgrade later!
      fData(x)$SNP10 = 0
      fData(x)$SNP10[ which(featureNames(x) %in% SNPs) ] = 1
    }
    betas(x)[ which(fData(x)$SNP10==1), ] = NA
    b = batch.id

    # ugly nasty kludge for DCC fiddling
    if( platform == 'HumanMethylation450' ) { 
      l3headers = c(rows='Composite Element REF', # {{{
                    Beta='Beta_value',
                    Pval='Detection_P_value',
                    Chrom='Chromosome',
                    Coord='Genomic_Coordinate') # }}}
    } else if( platform == 'HumanMethylation27') {
      l3headers = c(rows='Composite Element REF', # {{{
                    Beta='Beta_value',
                    Symbol='Gene_Symbol',
                    Chrom='Chromosome',
                    Coord='Genomic_Coordinate') # }}}
    }

    # reduce code duplication between serial & parallel 
    write.level3 <- function(s) { # {{{
      xs = which(sampleNames(x) == s)
      message(paste("Writing level 3 data for sample", which(subjects == s),
                    "of", length(subjects), "in TCGA batch", batchname))

      ## Level 3, 27k: masked betas, symbols, chr.hg18, site.hg18
      ## Level 3, 450k: masked betas, pvals, chr.hg18, site.hg18
      for( fvar in c('SYMBOL','CHR36','CPG36') ) { # for DCC use
        if( !(fvar %in% fvarLabels(x)) ) fData(x)[ , fvar ] = NA
      }
      lvl3data = data.frame( Beta=betas(x)[,xs], 
                             Pval=pvals(x)[,xs],            # dropped if 27k
                             Symbol=fData(x)[,'SYMBOL'],    # dropped if 450k
                             Chrom=fData(x)[,'CHR36'],      # must be hg18  
                             Coord=fData(x)[,'CPG36'] )     # must be hg18 
      lvl3data[ which(lvl3data$Pval > 0.05), 'Beta' ] <- NA # mask SNP/RPT too
      rownames(lvl3data) = featureNames(x)
      dump.file = paste(paste(diseasestub,platform,b,'lvl-3',s,'txt',sep='.'))
      headers1 = paste('Hybridization REF', s, s, sep="\t")
      headers2 = paste(l3headers, collapse="\t")
      cat(headers1, "\n", sep='', file=dump.file)
      cat(headers2, "\n", sep='', file=dump.file, append=TRUE)
      write.table(lvl3data[ , names(l3headers)[-1] ], # kludge
                  file=dump.file, quote=FALSE, append=TRUE,
                  row.names=TRUE, col.names=FALSE, sep="\t")
      return('done')
    } # }}}

    if(parallel) { 
      results <- mclapply(subjects, write.level3)
    } else { 
      results <- lapply(subjects, write.level3)
    }

    gc(,T)
    setwd(oldwd) 

  } # }}}

  message(paste(paste("Finished writing batch", batch.id, "level"), lvls))

} # }}}

## build aux, level 1, and mage-tab across all batches; levels 2 & 3 by batch
buildArchive<-function(map,version='0',base=NULL,platform='HumanMethylation450',                       lvls=c(1:3))
{ # {{{
  METHYLUMISET = FALSE
  if(is(map, 'MethyLumiSet')) { # {{{ METHYLUMISET = TRUE
    x = map
    map = pData(x)
    sampleNames(x) = map$TCGA.ID
    METHYLUMISET = TRUE 
  } # }}}
  stopifnot('TCGA.ID' %in% names(map))
  stopifnot('TCGA.BATCH' %in% names(map))
  if(!('BATCH.ID' %in% names(map))) { # {{{
    map$BATCH.ID = as.numeric(as.factor(map$TCGA.BATCH))
  } # }}}
  stopifnot('diseaseabr' %in% names(map))
  bs = seq_along(levels(as.factor(map$TCGA.BATCH)))
  message('Creating archive directories...')
  dirs = makeArchiveDirs(map, version=version, base=base, platform=platform)
  message('Writing level 2 and level 3 data...')
  if(METHYLUMISET==TRUE) {
    for(b in bs) {
      writeBatch(x, b, version=version, lvls=lvls)
    }
  } else {
    for(b in bs) {
      runBatchByID(map,b,base=base,platform=platform)
    }
  }
  message('Writing mage-tab IDF and SDRF files...')
  mageTab(map, version=version, base=base, platform=platform)
  message('Writing checksums for each directory...')
  justSign(map, base=base, platform=platform) 
  message('Validating the data archives...')
  validateArchive(map, base=base, platform=platform)
  message('If validation passed without errors, run packageAndSign, then SFTP.')
  # packageAndSign(map, base=base, platform=platform) 
} # }}}

## Adds MD5sums for all of the archive directories for a tumor. 
justSign <- function(map, base=NULL, platform='HumanMethylation450') 
{ # {{{
  oldwd = getwd()
  disease = unique(map$diseaseabr)
  stopifnot(length(disease) == 1)
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}

  dir.base = paste(base, 'tcga', disease, sep='/')
  setwd(dir.base)
  dir.patt = paste('^./jhu-usc.edu_',disease,'.',platform,sep='')
  dirs = grep(dir.patt, list.dirs(), value=TRUE)
  for(d in gsub('\\./','',dirs)) {
    setwd(d)
    if(file.exists('MANIFEST.txt')) file.remove('MANIFEST.txt')
    system('/usr/bin/env md5sum * > MANIFEST.txt')
    setwd(dir.base)
  }
  message(paste('Added manifests for each', disease,'directory.'))
  setwd(oldwd)
} # }}}

## Adds MD5sums and tarballs for all of the archive directories for a tumor. 
packageAndSign <- function(map, base=NULL, platform='HumanMethylation450') 
{ # {{{
  oldwd = getwd()
  disease = unique(map$diseaseabr)
  stopifnot(length(disease) == 1)
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}

  dir.base = paste(base, 'tcga', disease, sep='/')
  setwd(dir.base)
  dir.patt = paste('^./jhu-usc.edu_',disease,'.',platform,sep='')
  dirs = grep(dir.patt, list.dirs(), value=TRUE)
  for(d in gsub('\\./','',dirs)) {
    setwd(d)
    if(file.exists('MANIFEST.txt')) file.remove('MANIFEST.txt')
    system('/usr/bin/env md5sum * > MANIFEST.txt')
    system('/usr/bin/env md5sum -c MANIFEST.txt')  # check it!
    setwd(dir.base)
    d.tgz = paste(d, 'tar', 'gz', sep='.')
    system(paste('/usr/bin/env tar czf', d.tgz, d))
    system(paste('/usr/bin/env md5sum',d.tgz,'>',paste(d.tgz,'md5',sep='.')))
    system(paste('/usr/bin/env md5sum -c', paste(d.tgz,'md5',sep='.')))
  }
  message(paste('Added MD5sums and tarballs for', disease))
  setwd(oldwd)
} # }}}

validateArchive <- function(map,base=NULL,platform='HumanMethylation450',full=F)
{ # {{{
  oldwd=getwd()
  disease = unique(map$diseaseabr)
  stopifnot(length(disease) == 1)
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}
  message('Assuming $HOME/TCGA points to the validator directory...')
  archive = paste(base, 'tcga', disease, sep='/')
  setwd(archive)
  # so that we can get the validator to run in bypass mode...
  system('for i in `ls -d jhu* | grep -v gz`; do touch $i.tar.gz; done')
  if(full) { 
    validator.command = paste('./validate.sh', archive, '-centertype CGCC')
  } else { 
    validator.command = paste('./validate.sh', archive, # don't expand
                              '-centertype CGCC -bypass -noremote')
  }
  message(paste('Attempting to run', validator.command))
  setwd(paste( Sys.getenv('HOME'), 'TCGA', sep='/' ) )
  system(validator.command)
  setwd(oldwd) 
} # }}}
