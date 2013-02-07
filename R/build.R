## FIXME: the new world order
##
## Level 0 / AUX: processing logs, scripts, .SDF files, mappings
## Level 1: IDAT files ( two lines per sample in the SDRF )
## Level 2: background-corrected probe intensities ( M, U, Pval )
## Level 3: SNP10/pval-masked beta values with symbols, chrom.hg18, coord.hg18

## load and canonicalize
#loadMap <- function(disease, base=NULL, platform='HumanMethylation450') 
#{ # {{{ 
#  oldwd = getwd()
#  if(is.null(base)) { # {{{
#    message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')
#    if(platform == 'HumanMethylation27') {
#      base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # default
#    } else { 
#      base = paste( Sys.getenv('HOME'), 'meth450k', sep='/' ) # default
#    }
#  } # }}}
#  setwd( paste(base, 'mappings.aux', sep='/') )
#  map = read.csv( paste(disease, 'mappings', 'csv', sep='.'), 
#                  stringsAsFactors=F,
#                  row.names=1)
#  map = canonicalizeMapping(map)
#  map$BATCH.ID = as.character(as.numeric(as.factor(map$TCGA.BATCH)))
#  setwd(oldwd)
#  return(map)
#} # }}}

loadMap <- function(con, disease, platform="HumanMethylation450"){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	query <- paste("SELECT s.barcode, LOWER(s.uuid) AS uuid, s.samplename AS 'TCGA.ID', b.batch AS 'TCGA.BATCH', b.ordering AS 'BATCH.ID', h.name AS histology, t.name AS tissue, d.name AS diseaseabr, s.plate AS plate, s.well AS well FROM SAMPLE s INNER JOIN BATCH b ON b.id = s.batch INNER JOIN HISTOLOGY h ON h.id = s.histology INNER JOIN DISEASE d ON d.id = h.disease INNER JOIN TISSUE t ON t.id = h.tissue INNER JOIN STATUS st ON st.id = s.status INNER JOIN PLATFORM p ON p.id = s.platform INNER JOIN PROJECT pr WHERE d.name LIKE '%", disease, "%' AND st.id = 4 AND pr.id = 1 AND p.name = '", platform, "'", sep="")
	map <- dbGetQuery(con, query)
	cntl.index <- grep("-", map$diseaseabr, fixed=TRUE)
	tmp <- map[-cntl.index, ]
	if(length(cntl.index != 0)){
		cntl <- strsplit(map$diseaseabr[cntl.index], "-")
		cntl.disease.pos <- lapply(cntl, function(x){which(!is.na(match(x, disease)))})
		cntl <- strsplit(map$TCGA.BATCH[cntl.index], "-")
		cntl.batch <- list()
		cntl.ordering <- list()
		for(i in seq_along(cntl)){
			cntl.batch[[i]] <- cntl[[i]][cntl.disease.pos[[i]]]
			cntl.ordering[[i]] <- unique(tmp$BATCH.ID[which(tmp$TCGA.BATCH == cntl.batch[[i]])])
		}
		cntl.batch <- unlist(cntl.batch, use.names=F)
		cntl.ordering <- unlist(cntl.ordering, use.names=F)
		map$TCGA.BATCH[cntl.index] <- cntl.batch
		map$diseaseabr[cntl.index] <- disease
		#query <- paste("SELECT batch, ordering FROM BATCH WHERE batch IN ('", paste(map$TCGA.BATCH[cntl.index], collapse="','"), "')", sep="")
		#batch.cntl <- dbGetQuery(con, query)
		map$BATCH.ID[cntl.index] <- cntl.ordering
	}
	return(map)
}

## process a tumor's worth of Excel files from Dan
mapBatches <- function(xls.files, parallel=FALSE, link.raw=FALSE) 
{ # {{{
  if(parallel) {
    require(parallel)
    map = do.call(rbind, mclapply(xls.files, canonicalizeExcelFile, check=F))
  } else { 
    map = do.call(rbind, lapply(xls.files, canonicalizeExcelFile, check=F))
  }
  if(link.raw == TRUE) linkRawData(map, unlink.old.files=FALSE)
  map <- canonicalizeMapping(map)
  map$BATCH.ID = as.numeric(as.factor(map$TCGA.BATCH))
  return(map)
} # }}}

## Get MANIFEST file from TCGA portal
getManifest <- function(version='0', level='1', disease, batch, platform="HumanMethylation450"){
	target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
	domain <- "cgcc/jhu-usc.edu/humanmethylation450/methylation"
	center <- "jhu-usc.edu"
	dir <- paste("Level", level, sep="_")
	url <- paste(target,
		     tolower(disease),
		     domain,
		     paste(paste(center, disease, sep="_"), platform, dir, batch, version, 0, sep="."),
		     "MANIFEST.txt",
		     sep="/")
	cmd <- paste("curl -O", url, sep=" ")
	system(cmd)
}

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

## write level 1/2/3 for a batch -- aux/mage-tab are still done with a full map
writeBatch <- function(x,batch.id,old.version='0',version='0',base=NULL,parallel=F,lvls=c(1:3), revision=FALSE)
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
  if(!identical(sampleNames(x), x$TCGA.ID)) sampleNames(x) <- x$TCGA.ID
  dirs$raw = paste(base, 'raw', disease, sep='/')
  pkged = dirs$archive = paste(base, 'tcga', disease, paste("version", version, sep=""), sep='/')
  message('Assuming symlinks $HOME/meth27k and $HOME/meth450k both exist...')

  if(platform == 'HumanMethylation27') { # {{{
    require(IlluminaHumanMethylation27k.db)
    base = paste( Sys.getenv('HOME'), 'meth27k', sep='/' ) # }}}
  } else { # {{{
    require(IlluminaHumanMethylation450k.db)
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
    if(!revision){
      for(s in subjects) { # {{{
        xs = which(sampleNames(x) == s)
        message(paste("Copying IDAT files for sample", which(subjects == s),
                      "of", length(subjects), "in TCGA batch", batchname))
        stub = paste(x$barcode[xs], '_', sep='')
        for(i in paste(stub, c('Grn.idat','Red.idat'), sep='')) {
          file.link(paste(dirs$raw, i, sep='/'), paste('.', i, sep='/'))
        }
      } # }}}
    } else{
	message(paste("Retrieving Level 1 MANIFEST for Batch", batch.id, "from TCGA portal for revision", sep=" "))
        #getManifest(version=old.version, level='1', disease=disease, batch=batch.id, platform)
	#manifest <- paste(pkged,paste(diseasestub, platform, 'Level_1', batch.id,
        #                              version, '0', sep='.'), 'MANIFEST.txt', sep='/')
	#cmd <- paste('cp', manifest, '.', sep=' ')
	#system(cmd)
    }
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
    if(!revision){
      write.level2 <- function(s) { # {{{
        xs = which(sampleNames(x) == s)
        message(paste("Writing level 2 data for sample", which(subjects == s),
                      "of", length(subjects), "in TCGA batch", batchname))
        lvl2data = data.frame( M=methylated(x)[,xs], 
                               U=unmethylated(x)[,xs],
                               P=pvals(x)[,xs] )
        rownames(lvl2data) = featureNames(x)
        dump.file = paste(paste(diseasestub,platform,b,'lvl-2',s,'txt',sep='.'))
        headers1 = paste('Hybridization REF', s, s, s, sep="\t")
        headers2 = paste('Composite Element REF', 'Methylated_Intensity', 
                         'Unmethylated_Intensity', 'Detection_P_value', sep="\t")
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
  } else{
	message(paste("Retrieving Level 2 MANIFEST for Batch", batch.id, "from TCGA portal for revision", sep=" "))
        #getManifest(version=old.version, level='2', disease=disease, batch=batch.id, platform)
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
    b = batch.id

    # un-mask all the intensities (add pvals?!?)
    betas(x) <- methylated(x)/total.intensity(x)
    #if( !('SNP10' %in% fvarLabels(x)) ) {
    #  data(SNPs) # can upgrade later!
    #  fData(x)$SNP10 = 0
    #  fData(x)$SNP10[ which(featureNames(x) %in% SNPs) ] = 1
    #}
    if( !('mask' %in% fvarLabels(x))){
	    ifelse(platform == 'HumanMethylation27k', data(probesToMask.27k), data(probesToMask))
	    #data(probesToMask)
	    fData(x)$mask = 0
	    fData(x)$mask[ which(featureNames(x) %in% names(toMask))] = 1
    }
    #betas(x)[ which(fData(x)$SNP10==1), ] = NA
    betas(x)[ which(fData(x)$mask==1), ] = NA
    l3headers = c('Composite Element REF','Beta_value',
                  'Gene_Symbol','Chromosome','Genomic_Coordinate')
    #FIXME : Need to merge 27k annotations based on hg19 to level3.symbols.hg19
    additional.columns = l3headers[3:5] # as on the above line
    if(!all(additional.columns %in% fvarLabels(x))) { # {{{
      data(level3.symbols.hg19) # merged from 450k and 27k manifests
      missing.cols = setdiff(additional.columns, fvarLabels(x))
      fData(x) <- cbind(fData(x), level3.symbols.hg19[featureNames(x), missing.cols])
    } # }}}
  
    # reduce code duplication between serial & parallel 
    write.level3 <- function(s) { # {{{
      xs = which(sampleNames(x) == s)
      message(paste("Writing level 3 data for sample", which(subjects == s),
                    "of", length(subjects), "in TCGA batch", batchname))
      lvl3data = data.frame(Beta=betas(x)[,xs], 
                            Gene_Symbol=fData(x)[,'Gene_Symbol'],
                            Chromosome=fData(x)[,'Chromosome'],
                            Genomic_Coordinate=as.integer(fData(x)[,'Genomic_Coordinate']),
                            Pval=pvals(x)[,xs])
      rownames(lvl3data) = featureNames(x)
      lvl3data[ which(lvl3data$Pval > 0.05), 'Beta' ] <- NA
      lvl3data$Pval <- NULL # no longer included in Level 3 output
      lvl3data$Chromosome[which(is.na(lvl3data$Chromosome))] <- NA
      lvl3data$Genomic_Coordinate[which(is.na(lvl3data$Genomic_Coordinate))] <-0
      dump.file = paste(paste(diseasestub,platform,b,'lvl-3',s,'txt',sep='.'))
      headers1 = paste('Hybridization REF', s, s, s, s, sep="\t")
      headers2 = paste(l3headers, collapse="\t")
      cat(headers1, "\n", sep='', file=dump.file)
      cat(headers2, "\n", sep='', file=dump.file, append=TRUE)
      write.table(lvl3data,
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

  message(paste(paste("Wrote", disease, "batch", batch.id, "level"), lvls))

} # }}}

## build aux, level 1, and mage-tab across all batches; levels 2 & 3 by batch
buildArchive<-function(map, old.version='0', new.version='0', base=NULL,platform='HumanMethylation450', magetab.version=NULL, write.magetab=TRUE, lvls=c(1:3), revision=FALSE)
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
  ## FIXME: accomodate differing logic for LAML 450k vs. 27k
  #if('BATCH.ID' %in% names(map) && length(levels(as.factor(map$BATCH.ID)))==1){
  #  bs = unique(map$BATCH.ID)
  #} else {
  #  bs = seq_along(levels(as.factor(map$TCGA.BATCH)))
  #}
  bs = unique(map$BATCH.ID)
  bs = bs[order(bs)]
  if(is.null(magetab.version)){
    magetab.version = new.version
  }
  lvl = c("aux", "mage-tab")
  if(1 %in% lvls) lvl = c(lvl, "Level_1")
  if(2 %in% lvls) lvl = c(lvl, "Level_2")
  if(3 %in% lvls) lvl = c(lvl, "Level_3")
  #if(all(c(1:3) %in% lvls)) lvl = c("aux" , lvl)
  message('Creating archive directories...')
  dirs = makeArchiveDirs(map, version=new.version, base=base, magetab.version, lvls=lvl, platform=platform)
  writeTidyHistory(x, filepath=dirs$aux)
  message('Writing level 2 and level 3 data...')
  if(METHYLUMISET==TRUE) {
    for(b in bs) {
      writeBatch(x, b, old.version=old.version, version=new.version, lvls=lvls, revision=revision)
    }
  } else {
    for(b in bs) {
      runBatchByID(map,b,base=base,platform=platform)
    }
  }
  if(write.magetab){
    message('Writing mage-tab IDF and SDRF files...')
    mageTab(map, old.version=old.version, new.version=new.version, base=base, magetab.version=magetab.version, platform=platform, lvls=lvls, revision=revision)
  }
  #message('Packaging and signing each directory...')
  #packageAndSign(x, base=base, version=new.version, platform=platform) 
  #message('Validating the data archives...')
  #validateArchive(x, base=base, version=new.version, platform=platform)
  #message('If validation passed without errors, validate again on epigraph without -noremote, then SFTP.')
  # packageAndSign(map, base=base, platform=platform) 
} # }}}

## Adds MD5sums for all of the archive directories for a tumor. 
justSign <- function(map, base=NULL, version='0', platform='HumanMethylation450') 
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

  dir.base = paste(base, 'tcga', disease, paste("version", version, sep=""), sep='/')
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
packageAndSign <- function(map, base=NULL, version='0', platform='HumanMethylation450', revision=FALSE) 
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

  dir.base = paste(base, 'tcga', disease, paste("version", version, sep=""), sep='/')
  setwd(dir.base)
  dir.patt = paste('^./jhu-usc.edu_',disease,'.',platform,sep='')
  dirs = grep(dir.patt, list.dirs(), value=TRUE)
  for(d in gsub('\\./','',dirs)) {
    if(!(revision & grepl("Level_[12]", d))){
      setwd(d)
      if(file.exists('MANIFEST.txt')) file.remove('MANIFEST.txt')
      system('/usr/bin/env md5sum * > MANIFEST.txt')
      system('/usr/bin/env md5sum -c MANIFEST.txt')  # check it!
    }
    setwd(dir.base)
    d.tgz = paste(d, 'tar', 'gz', sep='.')
    system(paste('/usr/bin/env tar czf', d.tgz, d))
    system(paste('/usr/bin/env md5sum',d.tgz,'>',paste(d.tgz,'md5',sep='.')))
    system(paste('/usr/bin/env md5sum -c', paste(d.tgz,'md5',sep='.')))
  }
  message(paste('Added MD5sums and tarballs for', disease))
  setwd(oldwd)
} # }}}

## Does what it says on the tin -- tries to validate a newly built tumor archive
validateArchive <- function(map,base=NULL, version='0', platform='HumanMethylation450',full=F)
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
  archive = paste(base, 'tcga', disease, paste("version", version, sep=""), sep='/')
  setwd(archive)
  # so that we can get the validator to run in bypass mode...
  #system('for i in `ls -d jhu* | grep -v gz`; do touch $i.tar.gz; done')
  if(full) { 
    validator.command = paste('./validate.sh', archive, '-centertype CGCC')
  } else { 
    validator.command = paste('./validate.sh', archive, # don't expand
                              '-centertype CGCC -noremote')
  }
  message(paste('Attempting to run', validator.command))
  setwd(paste( Sys.getenv('HOME'), 'TCGA', sep='/' ) )
  system(validator.command)
  setwd(oldwd) 
} # }}}

## Wrapper function that makes building an archive basically a one liner
packageArchive <- function(map, disease=NULL, old.version='0', new.version='0', platform='HumanMethylation450',
			   magetab.version=NULL, write.magetab=TRUE, lvls=c(1:3), revision=FALSE)
{
	if(is.null(disease)){
		stop("Please provide the tumor name for which the archive is to be built")
	}
	
	cwd <- getwd()
	wd <- paste("~", disease, sep="/")
	raw <- paste("~", "meth450k", "raw", disease, sep="/")
	mset <- paste("~", "meth450k", "MethyLumiSets", sep="/")
	
	setwd(raw)
	
	message(paste("Reading idats for", disease, sep=" "))
		
	TUMOR <- methylumIDAT(map, parallel=T)

	message("Checking for Failed Samples")
	
	failed <- runSampleQC(TUMOR, filepath=wd)
	if(!is.null(failed)){
		TUMOR <- TUMOR[, -failed]
	}
	save(TUMOR, file=paste(mset, paste(disease, "raw", "rda", sep="."), sep="/"))
	gc()

	message("Performing Background Correction and Stripping MethyLumiSet of unneeded data")
	
	TUMOR <- stripMethyLumiSet(methylumi.bgcorr(TUMOR))
	gc()

	message("Performing Dye-Bias Equalization")
	
	TUMOR <- normalizeMethyLumiSet(TUMOR)
	save(TUMOR, file=paste(mset, paste(disease, "rda", sep="."), sep="/"))
	gc()

	message("Generating QC Probe Plot")
	
	pdf(paste(wd, paste(disease, "pdf", sep="."), sep="/"))
	qc.probe.plot(TUMOR)
	dev.off()
	gc()

	message("Generating Histogram of No. of Failed probes per sample")
	
	png(paste(wd, paste(disease, "sample", "summary", "png", sep="."), sep="/"))
	plotSampleSummary(TUMOR)
	dev.off()
	gc()

	message("Generating Density Plot of Cell Line Control Beta Values")
	
	controls <- which(TUMOR$histology %in% c("Cell Control Line", "Cytogenetically Normal", "Cell Line Control"))
	png(paste(wd, paste(disease, "cell", "line", "controls", "png", sep="."), sep="/"))
	plotDensities(TUMOR, controls = controls, label = "Cell Line Controls Beta")
	dev.off()
	gc()

	message("Writing out failure rate for Samples and Probes")
	
	writeSampleSummary(TUMOR, filepath = wd)
	writeProbeSummary(TUMOR, filepath = wd)
	gc()

	message("Building Level 1, 2, 3, aux and magetab Archives")
	
	buildArchive(TUMOR, base = base, old.version = old.version, new.version = new.version, platform = platform,
		     magetab.version = magetab.version, write.magetab = write.magetab, lvls = lvls, revision = revision)
	gc()

	message("Packaging archives into tarballs and generating MD5 sums")
	
	packageAndSign(TUMOR, base = base, version = new.version, platform = platform, revision = revision)

	setwd(cwd)
}
