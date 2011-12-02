getGlob <- function(barcodes) { # {{{
  glob.chip = substr(barcodes, 1, 10)
  glob.file = paste(barcodes, "_*.idat", sep='')
  glob = paste(paste('..', glob.chip, glob.file, sep='/'), collapse=' ')
  glob
} # }}}

map.IDs <- function(user=NULL, pass=NULL, ids=NULL, tissue=NULL, URL='http://webapp.epigenome.usc.edu/ECCP/samplemap.jsp') { # {{{

  require('RCurl')
  # sanity check: does it have 5 columns (no more, no less)?
  mapnames  <- c('beadchip','slot','ID','sex','tissue')

  if(is.null(user)) 
    user <- options('user')
  stopifnot(!is.null(user))
  if(is.null(pass)) 
    pass <- options('pass')
  stopifnot(!is.null(pass))

  beadmap <- strsplit(getURL(URL, userpwd=paste(user,pass,sep=':')),'\n')[[1]]
  beadmap <- beadmap[grep('\t',beadmap)]

  mapper <- sapply(beadmap, function(x) strsplit(x,'\t')[[1]])
  names(mapper) <- lapply(mapper, function(x) paste(x[[1]],x[[2]],sep='_'))
  mapper <- mapper[which(sapply(mapper, length)==length(mapnames))]
  mapper <- t(as.data.frame(mapper))

  colnames(mapper) <- mapnames
  rownames(mapper) <- gsub('^X','',rownames(mapper))
  mapper <- as.data.frame(mapper)

  mapper$slot <- as.ordered(mapper$slot)
  mapper$sex <- as.factor(mapper$sex)
  mapper$tissue <- as.factor(mapper$tissue)
  mapper$array <- rownames(mapper)

  if(!is.null(ids)) {
    mapper[ which(mapper$ID %in% ids), ]
  } else if(!is.null(tissue)) {
    mapper[ which(tolower(mapper$tissue) == tolower(tissue)), ]
  } else {
    mapper
  }

} # }}}

map.pheno <- function(user=NULL,pass=NULL,batch=NULL,platform=NULL, tumor=NULL, URL='http://webapp.epigenome.usc.edu/ECCP/arraymap.jsp') { # {{{

  jsondata <- map.json(user, pass, batch, platform, tumor, URL)
  browser() 
  mapper <- sapply(beadmap, function(x) strsplit(x,'\t')[[1]])

} # }}}

map.json <- function(user=NULL, pass=NULL, batch=NULL, tcgaID=NULL, disease=NULL, tissue=NULL, is.450k=TRUE, URL='http://webapp.epigenome.usc.edu/ECCP/arraymap.jsp') { # {{{

  if(is.null(user) && !is.null(options("user")[[1]])) user = options('user')
  if(is.null(pass) && !is.null(options("pass")[[1]])) pass = options('pass')
  stopifnot(!is.null(user) && !is.null(pass))
  userpwd = paste(user, pass, sep=':')
  require('RJSONIO')
  require('RCurl')
  jsonhttp = strsplit(getURL(URL, userpwd=userpwd), "\n")[[1]]
  jsondata = lapply(jsonhttp, function(x) {
    res = try(fromJSON(x))
    if(class(res) == 'try-error') return(NULL)
    else return(res)
  })
  parsed = chips.parse(jsondata, is.450k=is.450k)
  dat = parsed
  if( !is.null(tcgaID)) {
    if(length(tcgaID)>1) dat = dat[which(dat$tcgaID %in% tcgaID), ]
    else dat = dat[which(dat$tcgaID == tcgaID), ]
  }
  if( !is.null(disease)) {
    disease = tolower(disease)
    if(length(disease)>1) dat = dat[which(tolower(dat$disease) %in% disease),]
    else dat = dat[ which(tolower(dat$disease) == disease), ]
  }
  if( !is.null(batch)) {
    if(length(batch)>1) dat = dat[which(dat$batch %in% batch), ]
    else dat = dat[ which(dat$batch == batch), ]
  }
  if( !is.null(tissue)) {
    tissue = tolower(tissue)
    if(length(tissue)>1) dat = dat[which(tolower(dat$tissue) %in% tissue), ]
    else dat = dat[ which(tolower(dat$tissue) == tissue), ]
  }
  dat$barcode = rownames(dat)
  return(dat)
  
} # }}}

slots.27k.to.450k <- function(slots) { # {{{
  rs = paste('R0', 1:6, sep='')
  slotses = c(paste(rs, 'C01', sep=''), paste(rs, 'C02', sep=''))
  names(slotses) = toupper(letters[1:length(slotses)])
  return(slotses[toupper(slots)])
} # }}}

tcga.ID.from.name <- function(x) { # {{{
  res = strsplit(x, '-')[[1]]
  if(res[1] == 'TCGA') {
    return(as.character(res[[3]]))
  } else {
    return(NA)
  }
} # }}}

lanes.parse <- function(x, is.450k=TRUE) { # {{{
  nameses = names(x[['lane']][[1]])
  if(is.null(nameses)) return(NULL)
  dat = t(as.data.frame(lapply(x[['lane']], function(x) unlist(x))))
  colnames(dat) = nameses
  dat = as.data.frame(dat)
  dat$limsID = x[['flowcellProperties']][['limsID']]
  dat$tcgaID = sapply(levels(dat$name)[dat$name], tcga.ID.from.name)
  dat$slot27k = substr(dat$lane, 1, 1)
  dat$slot450k = slots.27k.to.450k(dat$slot27k)
  barcode = x[['flowcellProperties']][['serial']]
  if(is.450k) rownames(dat) = paste(barcode, dat$slot450k, sep='_')
  else rownames(dat) = paste(barcode, dat$slot27k, sep='_')
  return(dat)
} # }}}

chips.parse <- function(lst, is.450k=TRUE) { # {{{
  chs = lapply(lst, function(l) {
    res = try( lanes.parse(l, is.450k=is.450k) )
    if(class(res) == 'try-error') return(NULL)
    else return(res)
  })
  dat = chs[[1]]
  for(ch in 2:length(chs)) if(!is.null(chs[[ch]])) dat = rbind(dat, chs[[ch]])
  return(as.data.frame(dat))
} # }}}
