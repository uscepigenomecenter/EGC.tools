setClass("tcgaExtract",
	 representation(name="character",
			barcode = "character",
			protocol.label="Rle",
			name.label="character",
			label="character",
			term.source.mged="Rle",
			protocol.hybridize="Rle",
			name.hybridize="character",
			array.ref="Rle",
			term.source="Rle",
			protocol.image="Rle"))

setClass("tcgaLevel1",
	 representation(protocol.extract="Rle",
			name.scan="character",
			idat.file="character",
			name.archive.l1="Rle",
			data.type.l1="Rle",
			data.level.l1="Rle",
			include.l1="Rle"))

setClass("tcgaLevel2",
	 representation(protocol.norm="Rle",
			name.norm="character",
			data.matrix.file.l2="character",
			name.archive.l2="Rle",
			data.type.l2="Rle",
			data.level.l2="Rle",
			include.l2="Rle"))

setClass("tcgaLevel3",
	 representation(protocol.mask="Rle",
			name.mask="character",
			data.matrix.file.l3="character",
			name.archive.l3="Rle",
			data.type.l3="Rle",
			data.level.l3="Rle",
			include.l3="Rle"))


setClass("SDRF", contains=c("tcgaExtract","tcgaLevel1", "tcgaLevel2", "tcgaLevel3"),
	 representation(headers="character",
			sdrf.name="character",
			disease="character",
			platform="character",
			version="character"))

SDRF <- function(x, old.version='0', new.version='0', platform='HumanMethylation450'){
	if(is(x, 'MethyLumiSet')) {
		stopifnot('barcode' %in% varLabels(x))
		stopifnot('BATCH.ID' %in% varLabels(x))
		stopifnot('uuid' %in% varLabels(x))
		subjects = sampleNames(x) = x$TCGA.ID
		uuid = x$uuid
		platform = gsub('k$','',gsub('Illumina','',annotation(x)))
	} else {
		stopifnot('barcode' %in% names(x))
		stopifnot('BATCH.ID' %in% names(x))
		stopifnot('uuid' %in% names(x))
		subjects = rownames(x) = x$TCGA.ID
		uuid = x$uuid
	}

	diseases = levels(as.factor(x$diseaseabr))
	if(length(diseases)!=1) {
		stop('You have multiple disease abbreviations in this dataset.')
	} else {
		disease = diseases[1]
	}
	domain = 'jhu-usc.edu'
	prepreamble =  paste('jhu-usc.edu_',disease,'.',platform,sep='')
	preamble =  paste(prepreamble,'1',new.version,'0',sep='.')
	# <Domain>_<TumorType>.<Platform>.<ArchiveSerialIndex>.sdrf.txt
	sdrf.name = paste(preamble, 'sdrf','txt',sep='.') # all one big SDRF file
	array.ref = paste('Illumina.com','PhysicalArrayDesign',platform, sep=':')
	protocols = c(protocol1='labeling',
	              protocol2='hybridization',
	              protocol3='image_acquisition',
	              protocol4='feature_extraction',
	              protocol5='within_bioassay_data_set_function',
	              protocol6='within_bioassay_data_set_function')
	protocols = paste(domain, protocols,platform,'01',sep=':')
	termsources = c(termsource1='MGED Ontology',
	                termsource2='caArray')

	headers = c(extract='Extract Name',                         # uuid  {{{
	            barcode='Comment [TCGA Barcode]',               # TCGA ID
                    protocol1='Protocol REF',                       # labeling
                    labeled='Labeled Extract Name',                 # TCGA ID
                    label='Label',                                  # Cy3 or Cy5
                    termsource1='Term Source REF',                  # MGED Ontology
                    protocol2='Protocol REF',                       # hybridization
                    hybridization='Hybridization Name',             # barcode
                    arraydesign='Array Design REF',                 # array.ref
                    termsource2='Term Source REF',                  # caArray

                    ## Level 1 : IDAT files

                    protocol3='Protocol REF',                       # imaging 
                    protocol4='Protocol REF',                       # extraction
                    name1='Scan Name',                              # TCGA ID
                    datafile1='Array Data File',                    # _(Red|Grn).idat
                    archive1='Comment [TCGA Archive Name]',         # which batch
                    datatype1='Comment [TCGA Data Type]',           # Methylation
                    datalevel1='Comment [TCGA Data Level]',         # Level 1, duh
                    include1='Comment [TCGA Include for Analysis]', # yes

                    ## Level 2 : background-corrected M, U

                    protocol5='Protocol REF',                       # bg correction
                    name2='Normalization Name',                     # TCGA ID
                    datafile2='Derived Array Data Matrix File',     # .lvl-2.TCGA.txt
                    archive2='Comment [TCGA Archive Name]',         # which batch
                    datatype2='Comment [TCGA Data Type]',           # Methylation
                    datalevel2='Comment [TCGA Data Level]',         # Level 2, duh
                    include2='Comment [TCGA Include for Analysis]', # yes

                    ## Level 3 : masked betas, pvals/symbols, chromosome, coordinate

                    protocol6='Protocol REF',                       # bg correction
                    name3='Normalization Name',                     # TCGA ID
                    datafile3='Derived Array Data Matrix File',     # .lvl-3.TCGA.txt
                    archive3='Comment [TCGA Archive Name]',         # which batch
                    datatype3='Comment [TCGA Data Type]',           # Methylation
                    datalevel3='Comment [TCGA Data Level]',         # Level 3, duh
                    include3='Comment [TCGA Include for Analysis]') # yes }}}

	name = rep(uuid, each=2)
	barcode = rep(subjects, each=2)
	protocol.label = Rle(paste(domain,'labeling',platform,'01',sep=':'), length(subjects) * 2)
	#name.label = name
	label = rep(c('Cy3', 'Cy5'), length(subjects))
	term.source.mged = Rle("MGED Ontology", length(subjects) * 2)
	protocol.hybridize = Rle(paste(domain,'hybridization',platform,'01',sep=':'), length(subjects) * 2)
	#name.hybridize = name
	array.ref = Rle(paste('Illumina.com','PhysicalArrayDesign',platform,sep=':'), length(subjects) * 2)
	term.source = Rle('caArray', length(subjects) * 2)
	protocol.image = Rle(paste(domain,'image_acquisition',platform,'01',sep=':'), length(subjects) * 2)

	extract <- new("tcgaExtract",
		       name=name,
		       barcode=barcode,
		       protocol.label=protocol.label,
		       name.label=name,
		       label=label,
		       term.source.mged=term.source.mged,
		       protocol.hybridize=protocol.hybridize,
		       name.hybridize=name,
		       array.ref=array.ref,
		       term.source=term.source,
		       protocol.image=protocol.image)

	protocol.extract = Rle(paste(domain,'feature_extraction',platform,'01',sep=':'), length(subjects) * 2)
	#name.scan = name
	barcode = rep(x$barcode, each=2)
	channel = rep(c('Grn.idat', 'Red.idat'), length(subjects))
	idat.file = paste(barcode, channel, sep='_')
	batch = rep(x$BATCH.ID, each=2)
	name.archive.l1 = paste(prepreamble,'Level_1',batch,new.version,'0',sep='.')
	name.archive.l1 = as(name.archive.l1, "Rle")
	data.type = Rle('DNA Methylation', length(subjects) * 2)
	data.level.l1 = Rle('Level 1', length(subjects) * 2)
	include = Rle('yes', length(subjects) * 2)

	level1 <- new("tcgaLevel1",
		      protocol.extract=protocol.extract,
		      name.scan=name,
		      idat.file=idat.file,
		      name.archive.l1=name.archive.l1,
		      data.type.l1=data.type,
		      data.level.l1=data.level.l1,
		      include.l1=include)

	protocol.norm = Rle(paste(domain,'within_bioassay_data_set_function',platform,'01',sep=':'), length(subjects) * 2)
	#name.norm = name
	data.matrix.file.l2 = paste(prepreamble,batch,'lvl-2',name,'txt',sep='.')
	name.archive.l2 = paste(prepreamble,'Level_2',batch,new.version,'0',sep='.')
	name.archive.l2 = as(name.archive.l2, "Rle")
	data.level.l2 = Rle('Level 2', length(subjects) * 2)

	level2 <- new("tcgaLevel2",
		      protocol.norm=protocol.norm,
		      name.norm=name,
		      data.matrix.file.l2=data.matrix.file.l2,
		      name.archive.l2=name.archive.l2,
		      data.type.l2=data.type,
		      data.level.l2=data.level.l2,
		      include.l2=include)

	protocol.mask = Rle(paste(domain,'within_bioassay_data_set_function',platform,'01',sep=':'), length(subjects) * 2)
	#name.mask = name
	data.matrix.file.l3 = paste(prepreamble,batch,'lvl-3',name,'txt',sep='.')
	name.archive.l3 = paste(prepreamble,'Level_3',batch,new.version,'0',sep='.')
	name.archive.l3 = as(name.archive.l3, "Rle")
	data.level.l3 = Rle('Level 3', length(subjects) * 2)

	level3 <- new("tcgaLevel3",
		      protocol.mask=protocol.mask,
		      name.mask=name,
		      data.matrix.file.l3=data.matrix.file.l3,
		      name.archive.l3=name.archive.l3,
		      data.type.l3=data.type,
		      data.level.l3=data.level.l3,
		      include.l3=include)

	sdrf <- new("SDRF",
		    headers=headers,
		    sdrf.name=sdrf.name,
		    disease=disease,
		    platform=platform,
		    version=new.version,
		    extract,
		    level1,
		    level2,
		    level3)
}

setMethod("show", signature(object="SDRF"),
	  function(object){
		  cat(class(object), "\n")
		  n <- length(object@name) / 2
		  cat("No. of samples :", n, "\n", sep=" ")
		  cat("Tumor type :", object@disease, "\n", sep=" ")
		  cat("Platform :", object@platform, "\n", sep=" ")
		  cat("Data Type :", "DNA Methylation", "\n", sep=" ")
		  cat("Data Levels :", "1,2,3", "\n", sep=" ")
		  cat("MageTab version :", object@version, "\n", sep=" ")
	  })

setGeneric("getExtract", function(object) standardGeneric("getExtract"))
setGeneric("getLevel1", function(object) standardGeneric("getLevel1"))
setGeneric("getLevel2", function(object) standardGeneric("getLevel2"))
setGeneric("getLevel3", function(object) standardGeneric("getLevel3"))

setMethod("getExtract", signature(object="SDRF"),
	  function(object){
		  slots <- slotNames("tcgaExtract")
		  #v <- paste("object", slots, sep="@")
		  extract <- data.frame(object@name,
					as.character(object@protocol.label),
					object@name.label,
					object@label,
					as.character(object@term.source.mged),
					as.character(object@protocol.hybridize),
					object@name.hybridize,
					as.character(object@array.ref),
					as.character(object@term.source),
					as.character(object@protocol.image))
		  colnames(extract) <- slots
		  return(extract)
	  })

setMethod("getLevel1", signature(object="SDRF"),
	  function(object){
		  slots <- slotNames("tcgaLevel1")
		  #v <- paste("object", slots, sep="@")
		  level1 <- data.frame(as.character(object@protocol.extract),
				       object@name.scan,
				       object@idat.file,
				       as.character(object@name.archive.l1),
				       as.character(object@data.type.l1),
				       as.character(object@data.level.l1),
				       as.character(object@include.l1))
		  colnames(level1) <- slots
		  return(level1)
	  })

setMethod("getLevel2", signature(object="SDRF"),
	  function(object){
		  slots <- slotNames("tcgaLevel2")
		  #v <- paste("object", slots, sep="@")
		  level2 <- data.frame(as.character(object@protocol.norm),
				       object@name.norm,
				       object@data.matrix.file.l2,
				       as.character(object@name.archive.l2),
				       as.character(object@data.type.l2),
				       as.character(object@data.level.l2),
				       as.character(object@include.l2))
		  colnames(level2) <- slots
		  return(level2)
	  })

setMethod("getLevel3", signature(object="SDRF"),
	  function(object){
		  slots <- slotNames("tcgaLevel3")
		  #v <- paste("object", slots, sep="@")
		  level3 <- data.frame(as.character(object@protocol.mask),
				       object@name.mask,
				       object@data.matrix.file.l3,
				       as.character(object@name.archive.l3),
				       as.character(object@data.type.l3),
				       as.character(object@data.level.l3),
				       as.character(object@include.l3))
		  colnames(level3) <- slots
		  return(level3)
	  })
