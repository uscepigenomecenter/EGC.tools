#Based on the capwords function demonstrated as an example in the help manual for tolower
capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s,1,1)),
                       {s <- substring(s,2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }

# Utility function to aggreagate and clean up Dan's mappings before adding them to DANEUS
processMap <- function(base=NULL, disease, platform='HumanMethylation450') {
	if(is.null(base)){
		message("Assuming there is a directory called DaneusMappings in ~/Dropbox which contains Dan's mapping files")
		base <- "~/Dropbox/DaneusMappings"
	}
	map.dir <- ifelse(grepl('450', platform), file.path(base, disease), file.path(paste(base, "27k", sep="/"), disease))
	oldwd <- getwd()
	setwd(map.dir)
	message("Canonicalizing Dan's mapping files assuming they have been standardized to expected input")
	map <- mapBatches(list.files(pattern='.xls'))
	setwd(base)
	cat(map$TCGA.ID, file="barcodes.txt", sep="\n")
	message("Grabbing UUIDs from the Biospecimen Metadata Browser and expecting getuuid.pl script to be present in base directory")
	system("perl getuuid.pl barcodes.txt uuid.txt", ignore.stdout = TRUE)
	setwd(map.dir)
	uuid <- read.delim(file.path(base, "uuid.txt"), stringsAsFactors=F, header=F)[,2]
	if(any(uuid %in% c("", " "))){
		stop("There are missing UUIDs. Please check and remove the samples with missing UUIDs before proceeding")
	}
	map$uuid <- uuid
	map <- map[,-8]
	source <- lapply(strsplit(map$source, split=" "), function(x){x[c(2,4)]})
	plate <- unlist(lapply(source, function(x){x[1]}), use.names=F)
	well <- unlist(lapply(source, function(x){x[2]}), use.names=F)
	map$plate <- plate
	map$well <- well
	map$source <- NULL
	map$status <- "Infinium"
	hist <- map$histology
	# Replace multiple spaces with one space
	hist <- gsub("\\s+", " ", hist, perl=T)
	# Capitalise first letter of every word
	hist <- sapply(hist, capwords, strict=T, USE.NAMES=F)
	tissue <- map$tissue
	# Remove any spaces that might be at the end
	tissue <- gsub("\\s$", "", tissue, perl=T)
	# Rename the tissue type for all Normals to Matched Normals
	tissue[which(hist %in% c("Normal Tissue", "Normal"))] <- "Matched Normal"
	tissue[which(hist %in% c("Cytogenetically Normal", "Cell Control Line", "Control Cell Line"))] <- "Cell Line Control"
	tissue[which(tissue %in% "Rectal")] <- "Rectum"
	hist[which(hist == "Normal")] <- "Normal Tissue"
	hist[which(hist %in% c("Cytogenetically Normal", "Cell Control Line", "Control Cell Line"))] <- "Cell Line Control"
	hist <- gsub("Rectum", "Rectal", hist)
	hist <- gsub("Gbm", "GBM", hist)
	hist <- gsub("Iii", "III", hist, fixed=T)
	hist <- gsub("Ii", "II", hist, fixed=T)
	hist <- gsub("aml", "AML", hist, fixed=T)
	hist <- gsub("Nos", "NOS", hist, fixed=T)
	hist <- gsub("nos", "NOS", hist, fixed=T)
	hist <- gsub("lms", "LMS", hist, fixed=T)
	hist <- gsub("ups", "UPS", hist, fixed=T)
	map$histology <- hist
	map$tissue <- tissue
	rm(hist, tissue, plate, well, uuid, source)
	message("Need to add shipping dates manually, will fix it to automatically grab shipping dates once DCC 2.0 goes live")
	setwd(oldwd)
	return(map)
}

# Bunch of sanity checks to run on a mapping file before inserting the mappings into DANEUS
runSanityChecks <- function(map){
	message("Listing the  unique disease types in the mapping \n")
	print(unique(map$diseaseabr))
	message("Listing the unique tissue types in the mapping \n")
	print(unique(map$tissue))
	message("Listing unique histologies. Cross-check with know histologies in DANEUS to avoid spelling errors and duplication \n")
	print(unique(map$histology))
	message("Listing the batches and their ordering. Check carefully \n")
	print(unique(map$TCGA.BATCH))
	print(unique(map$BATCH.ID))
	if(any(duplicated(map$barcode))){
		stop("Duplicate samples present. Please check mapping and remove any duplicates \n")
	}
	if(any(duplicated(map$TCGA.ID))){
		stop("Duplicate samples present. Please check mapping and remove any duplicates \n")
	}
	if(any(is.na(map$TCGA.BATCH))){
		stop("There are missing batches in the mapping \n")
	}
}
	
