## Adds an Investigation Design Format file for all the batches in an archive.
buildIDF <- function(x, version='0', platform='HumanMethylation450') 
{ # {{{
  if(is(x, 'MethyLumiSet')) {
    platform = gsub('k$','',gsub('Illumina','',annotation(x)))
  }
  disease = unique(x$diseaseabr)
  stopifnot(length(disease)==1)
  filenm=paste('jhu-usc.edu_',disease,'.',platform,'.1.',version,'.0.idf.txt',
               sep='')
  sdrfnm=gsub('idf','sdrf',filenm)

  invtitle = paste('TCGA Analysis of DNA Methylation for', disease, 'using Illumina Infinium', platform, 'platform')
  pub.date = paste('Public Release Date', date(), sep="\t")
  protocoln = c('labeling','hybridization','scan')
  protocols = c('labeling','hybridization','image_acquisition')
  protocols = paste('jhu-usc.edu',protocols,platform,'01',sep=':')
  names(protocols) = protocoln

  ## FIXME: break this up properly and/or build it up properly in a loop
  bogosity = c(   
    paste(c('Investigation Title',
            invtitle),
            collapse="\t"),
    '', # chunk separator
    paste(c('Experimental Design',
            'disease_state_design'),
            collapse="\t"),
    paste(c('Experimental Design Term Source REF',
            'MGED Ontology'),
            collapse="\t"),
    paste(c('Experimental Factor Type',
            'disease'),
            collapse="\t"),
    paste(c('Experimental Factor Type Term Source REF',
            'MGED Ontology'),
            collapse="\t"),
    '', # chunk separator
    paste(c('Person Last Name',
            'Laird'),
            collapse="\t"),
    paste(c('Person First Name',
            'Peter'),
            collapse="\t"),
    paste(c('Person Mid Initials',
            'W'),
            collapse="\t"),
    paste(c('Person Email',
            'plaird@usc.edu'),
            collapse="\t"),
    paste(c('Person Phone',
            '323.442.7890'),
            collapse="\t"),
    paste(c('Person Address',
      'USC Epigenome Center, University of Southern California, CA 90033, USA'),            collapse="\t"),
    paste(c('Person Affiliation',
            'University of Southern California'),
            collapse="\t"),
    paste(c('Person Roles',
            'submitter'),
            collapse="\t"),
    '', # chunk separator
    paste(c('Quality Control Types',
            'real_time_PCR_quality_control'),
            collapse="\t"),
    paste(c('Quality Control Types Term Source REF',
            'MGED Ontology'),
            collapse="\t"),
    paste(c('Replicate Type',
            'bioassay_replicate_reduction'),
            collapse="\t"),
    paste('Replicate Type Term Source REF','MGED Ontology',sep="\t"),
    # paste('Date of Experiment ', date(), sep="\t"),
    paste('Public Release Date', date(), sep="\t"),
    paste('Protocol Name', paste(protocols, collapse="\t"), sep="\t"),
    paste('Protocol Type', paste(names(protocols), collapse="\t"), sep="\t"),
    paste(c('Protocol Term Source REF',
            'MGED Ontology',
            'MGED Ontology',
            'MGED Ontology'), 
            collapse="\t"),
    #
    # FIXME: couple this to the SDRF 
    #           
    # protocols = c(protocol1='labeling',
    #               protocol2='hybridization',
    #               protocol3='image_acquisition',
    #               protocol4='feature_extraction',
    #               protocol5='within_bioassay_data_set_function',
    #               protocol6='within_bioassay_data_set_function') 
    #           
    paste(c('Protocol Description',
            paste('Illumina Infinium', platform, 'Labeled Extract'),
            paste('Illumina Infinium', platform, 'Hybridization Protocol'),
            paste('Illumina Infinium', platform, 'Scan Protocol')),
            collapse="\t"),
    'Protocol Parameters',
    '', # chunk separator
    paste(c('SDRF Files', 
            sdrfnm), 
            collapse="\t"),
    paste(c('Term Source Name',
            'MGED Ontology',
            'caArray'),
            collapse="\t"),
    paste(c('Term Source File',
            'http://mged.sourceforge.net/ontologies/MGEDontology.php',
            'http://caarraydb.nci.nih.gov/'),
            collapse="\t"),
    paste(c('Term Source Version',
            '1.3.1.1',
            '2007-01'),
            collapse="\t")
  )

  cat(bogosity, sep="\n", file=filenm)

} # }}}

## Adds Sample and Data Relationship Format file for all batches in an archive.
buildSDRF <- function(x, version='0',platform='HumanMethylation450',lvls=c(1:3))
{ # {{{

  if(is(x, 'MethyLumiSet')) { # {{{
    stopifnot('barcode' %in% varLabels(x))
    stopifnot('BATCH.ID' %in% varLabels(x))
    subjects = sampleNames(x) = x$TCGA.ID    
    platform = gsub('k$','',gsub('Illumina','',annotation(x))) # }}}
  } else { # {{{
    stopifnot('barcode' %in% names(x))
    stopifnot('BATCH.ID' %in% names(x))
    subjects = rownames(x) = x$TCGA.ID    
  } # }}}

  diseases = levels(as.factor(x$diseaseabr))
  if(length(diseases)!=1) {
    stop('You have multiple disease abbreviations in this dataset.')
  } else {
    disease = diseases[1]
  }
  prepreamble =  paste('jhu-usc.edu_',disease,'.',platform,sep='')
  preamble =  paste(prepreamble,'1',version,'0',sep='.')
  # <Domain>_<TumorType>.<Platform>.<ArchiveSerialIndex>.sdrf.txt
  sdrf.name = paste(preamble, 'sdrf','txt',sep='.') # all one big SDRF file
  array.ref = paste('Illumina.com','PhysicalArrayDesign',platform, sep=':')
  protocols = c(protocol1='labeling', # {{{
                protocol2='hybridization',
                protocol3='image_acquisition',
                protocol4='feature_extraction',
                protocol5='within_bioassay_data_set_function',
                protocol6='within_bioassay_data_set_function') # }}}
  protocols = paste('jhu-usc.edu',protocols,platform,'01',sep=':')
  termsources = c(termsource1='MGED Ontology',
                  termsource2='caArray')

  ## FIXME: there's got to be a better way of handling all these fields
  headers = c(extract='Extract Name',                         # TCGA ID  {{{
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

  elements= c(extract='TCGA.ID.HERE', # {{{
              protocol1=protocols[1],
              labeled='TCGA.ID.HERE',
              label='Cy3 or Cy5',
              termsource1=termsources[1],
              protocol2=protocols[2],
              hybridization='barcode.here',
              arraydesign=array.ref,
              termsource2=termsources[2],

              # Level 1
              protocol3=protocols[3],
              protocol4=protocols[4],
              name1='TCGA.ID.HERE',
              datafile1='barcode_Red.or.Grn.idat',
              archive1='TCGA.BATCH.ARCHIVENAME',
              datatype1='DNA Methylation',
              datalevel1='Level 1',
              include1='yes',

              # level 2
              protocol5=protocols[5],
              name2='TCGA.ID.HERE',
              datafile2='jhu-usc.edu_blablah-lvl2.TCGA.txt',
              archive2='TCGA.BATCH.ARCHIVENAME',
              datatype2='DNA Methylation',
              datalevel2='Level 2',
              include2='yes',

              # level 3
              protocol6=protocols[6],
              name3='TCGA.ID.HERE',
              datafile3='jhu-usc.edu_blablah-lvl3.TCGA.txt',
              archive3='TCGA.BATCH.ARCHIVENAME',
              datatype3='DNA Methylation',
              datalevel3='Level 3',
              include3='yes') # }}}

  cat(paste(headers, collapse="\t"), "\n", sep='', file=sdrf.name)
  channels = c(Cy3='Grn',Cy5='Red')
  for(s in subjects) {  # {{{ add two lines in the SDRF for each subject!
    elements['extract'] = s
    elements['labeled'] = s
    elements['name1'] = s
    elements['name2'] = s
    elements['name3'] = s
    xs = which(x$TCGA.ID == s)
    b = x$BATCH.ID[ xs ]
    chip = x$barcode[ xs ]

    # not likely to use these, but just in case...
    if( !('1' %in% lvls) ) elements['include1'] = 'no'
    if( !('2' %in% lvls) ) elements['include2'] = 'no'
    if( !('3' %in% lvls) ) elements['include3'] = 'no'

    # {{{ level 1: IDAT files (2 lines per sample)
    # datafile1 will be two entries: one for Cy5 (red) and one for Cy3 (green)
    elements['archive1'] = paste(prepreamble,'Level_1',b,version,'0', sep='.') 
    # }}}

    # {{{ level 2: background-corrected M and U intensities
    elements['datafile2'] = paste(prepreamble,b,'lvl-2',s,'txt',sep='.')
    message(paste('Subject',s,'Level 2 data file is ', elements['datafile2']))
    elements['archive2'] = paste(prepreamble,'Level_2',b,version,'0', sep='.')
    message(paste('Subject',s,'Level 2 data file is in', elements['archive2']))
    # }}}

    # {{{ level 3: SNP10-masked beta value with ECDF pvalue
    elements['datafile3'] = paste(prepreamble,b,'lvl-3',s,'txt', sep='.')
    message(paste('Subject',s,'Level 3 data file is ', elements['datafile3']))
    elements['archive3'] = paste(prepreamble,'Level_3',b,version,'0', sep='.')
    message(paste('Subject',s,'Level 3 data file is in', elements['archive3']))
    # }}}

    for(ch in names(channels)) { # {{{ each subject has TWO SDRF lines
      elements['label'] = ch
      elements['hybridization'] = s # otherwise 
      elements['datafile1'] = paste(chip,'_',channels[[ch]],'.idat', sep='')
      cat(paste(elements,collapse="\t"),"\n",sep='',file=sdrf.name,append=T)
    } # }}}

  } # }}} 

} # }}}

mageTab <- function(map, version='0', base=NULL, platform='HumanMethylation450')
{ # {{{
  if(is(map, 'MethyLumiSet')) { # {{{
    x <- map
    platform = gsub('^Illumina', '', gsub('k$','', annotation(x)))
    map <- pData(x)
  } # }}}
  disease = unique(map$diseaseabr)
  stopifnot(length(disease) == 1)
  archive.dir = paste(Sys.getenv('HOME'), 'meth450k', 'tcga', disease, 
                      paste(paste( 'jhu-usc.edu', disease, sep='_'), platform,
                            'mage-tab','1',version,'0',sep='.'),sep='/')
  if(!('BATCH.ID' %in% names(map))) { # {{{
    stopifnot('TCGA.BATCH' %in% names(map))
    map$BATCH.ID = as.numeric(as.factor(map$TCGA.BATCH))
  } # }}}
  oldwd = getwd()
  setwd(archive.dir)
  addDescription(map)
  buildIDF(map, version=version)
  buildSDRF(map, version=version)
  setwd(oldwd)
} # }}} 

addDescription <- function(x, platform='HumanMethylation450') {
  #{{{
  if(is(x, 'MethyLumiSet')) {
    platform = gsub('k$','',gsub('Illumina','',annotation(x)))
  }
  disease = unique(x$diseaseabr)
  stopifnot(length(disease) == 1)
  boilerplate = paste('This data archive contains the Cancer Genome Atlas (TCGA) analysis of DNA methylation profiling using the IIllumina Infinium',platform,'platform. The Infinium platform analyzes up to 482,421 CpG dinucleotides and 3091 CpH trinucleotides, spanning gene-associated elements as well as intergenic regions. DNA samples were received, bisulfite converted and cytosine methylation was evaluated using IIllumina Infinium',platform,'microarrays.')
  sample.desc = paste('This archive contains Infinium DNA methylation data for',
                      disease, 'samples.')
  level.desc = 'Data levels and the files contained in each data level package are as follows:'
  level.desc = paste(level.desc,'',
'AUX: Auxilary directory containing .sdf files, mappings, and processing logs.',
'LEVEL 1: Level 1 data contain raw IDAT files (two per sample) as produced by the iScan system and as mapped by the SDRF and also in the disease mapping file.',
'LEVEL 2: Level 2 data contain background-corrected methylated (M) and unmethylated (U) summary intensities as extracted by the methylumi package.  Background correction is performed using a normal-exponential convolution.',
'LEVEL 3: Derived summary measures (beta values: M/(M+U) for each locus) with non-detection probabilities (P-values) computed using the ECDF of the negative control probes in each color channel.  Probes annotated as having a SNP within 10bp of the interrogated site are masked as NA across all samples in a batch, and probes with a non-detection probability greater than 0.05 are also masked.',
  sep="\n")
  cat(boilerplate, sample.desc, level.desc, file='DESCRIPTION.txt', sep="\n")
} # }}}
