capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s,1,1)),
                       {s <- substring(s,2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }

cleanMap <- function(path=".", disease) {
	map.dir <- file.path(path, disease)
	oldwd <- getwd()
	setwd(map.dir)
	map <- mapBatches(list.files(pattern='.xls'))
	cat(map$TCGA.ID, file="barcodes.txt", sep="\n")
	system("perl getuuid.pl barcodes.txt uuid.txt")
}
	
