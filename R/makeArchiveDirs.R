## FIXME: make this work for HumanMethylation27 as well (add a 'platform' arg)
##
makeArchiveDirs <- function(map, version='0', base=NULL, lvls=c('aux','Level_1','Level_2','Level_3','mage-tab'), platform='HumanMethylation450') 
{ # {{{

  METHYLUMISET = FALSE
  if( is(map, 'MethyLumiSet') ) { # {{{
    x = map
    map = pData(x)
    METHYLUMISET = TRUE
    platform = gsub('k$','',gsub('Illumina','',annotation(x))) 
    stopifnot('TCGA.ID' %in% varLabels(x))
    sampleNames(x) <- x$TCGA.ID
    message("Warning: this function assumes you have already removed controls!")
  } # }}}
  stopifnot( is(map, 'data.frame') )
  if(any(map$histology == 'Cytogenetically Normal')) { # {{{ stop & exit
    stop('You have cell line controls in your data.  Remove them.')
  } # }}}
  dirs = list()
  if(is.null(base)) { # {{{
    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
    if(platform == 'HumanMethylation27') {
      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
    } else { 
      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
    }
  } # }}}
  disease = unique(map$diseaseabr)
  dirs$disease = paste(base, 'raw', disease, sep='/')
  dirs$mappings = paste(base, 'mappings.aux', sep='/')
  diseasestub = paste('jhu-usc.edu', disease, sep='_')
  dirs$archive = pkged = paste(base, 'tcga', disease, sep='/')
  setwd(dirs$disease)
  oldwd = getwd()
  names(lvls) = lvls
  batches = levels(as.factor(map$TCGA.BATCH))

  if('aux' %in% lvls) { # {{{
    dirs$aux = paste(pkged, 
                     paste(diseasestub, 'HumanMethylation450', 
                           'aux', '1', version, '0', sep='.'),  sep='/')
    oldwd = getwd()
    dir.create(dirs$aux)
    setwd(dirs$aux)
    diseasemap = paste(disease, 'mappings', 'csv', sep='.')
    if(METHYLUMISET) writeTidyHistory(x, dirs$aux)
    write.csv(map, file=diseasemap)
    setwd(oldwd)
    system(paste('touch ', dirs$aux, '.tar.gz', sep=''))
    lvls = lvls[-which(names(lvls)=='aux')]
  } # }}}
  if('Level_1' %in% lvls) { # {{{
    dirs$level_1 = c()
    for( i in 1:length(batches) ) {
      batchnum = i
      batchname = batches[i]
      batchmap = map[ which(map$TCGA.BATCH == batchname), ]
      oldwd = getwd()
      dirs$level_1[i] = paste(pkged, 
                              paste(diseasestub, 'HumanMethylation450', 
                                    'Level_1',i,version,'0', sep='.'), sep='/')
      dir.create(dirs$level_1[i])
      # setwd(dirs$level_1[i])
      # for( beadchip in unique(substr(batchmap$barcode, 1, 10))) {
      #   which.samples = which(substr(batchmap$barcode, 1, 10) == beadchip)
      #   cmd=paste('ln -f ',dirs$disease,'/',beadchip,'.sdf ',dirs$aux,sep='')
      #   system(command=cmd)
      #   for( chip in grep(beadchip, batchmap$barcode, value=T) ) {
      #     cmd = paste('ln -f ', dirs$disease, '/', chip, '_*.idat .', sep='')
      #     system(command=cmd)
      #   }
      # }
      # setwd(oldwd)
      system(paste('touch ', dirs$level_1[i], '.tar.gz', sep=''))
    }
  } # }}}
  if('Level_2' %in% lvls) { # {{{
    dirs$level_2 = c()
    for( i in 1:length(batches) ) {
      batchnum = i
      batchname = batches[i]
      batchmap = map[ which(map$TCGA.BATCH == batchname), ]
      oldwd = getwd()
      dirs$level_2[i] = paste(pkged, 
                              paste(diseasestub, 'HumanMethylation450', 
                                    'Level_2',i,version,'0', sep='.'), sep='/')
      dir.create(dirs$level_2[i])
      system(paste('touch ', dirs$level_2[i], '.tar.gz', sep=''))
    }
    # if(METHYLUMISET) level2(x, version=version)
  } # }}}
  if('Level_3' %in% lvls) { # {{{
    dirs$level_3 = c()
    for( i in 1:length(batches) ) {
      batchnum = i
      batchname = batches[i]
      batchmap = map[ which(map$TCGA.BATCH == batchname), ]
      oldwd = getwd()
      dirs$level_3[i] = paste(pkged, 
                              paste(diseasestub, 'HumanMethylation450', 
                                    'Level_3',i,version,'0', sep='.'), sep='/')
      dir.create(dirs$level_3[i])
      system(paste('touch ', dirs$level_3[i], '.tar.gz', sep=''))
    }
    if(METHYLUMISET) level3(x, version=version)
  } # }}}
  if('mage-tab' %in% lvls) { # {{{
    dirs$magetab = paste(pkged, 
                         paste(diseasestub, 'HumanMethylation450', 
                               'mage-tab', '1', version, '0', sep='.'), sep='/')
    dir.create(dirs$magetab)
    system(paste('touch ', dirs$magetab, '.tar.gz', sep=''))
    lvls = lvls[-which(names(lvls)=='mage-tab')]
  } # }}}

  return(dirs)

} # }}}

