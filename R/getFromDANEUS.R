#A bunch of get methods to retrieve different information from DANEUS

getDisease <- function(con=NULL){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	disease <- dbGetQuery(con, "SELECT * FROM DISEASE")
	return(disease)
}

insertDisease <- function(con=NULL, disease){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	query <- paste("INSERT INTO DISEASE (name) VALUES ('", disease, "')", sep="")
	t <- tryCatch(dbSendQuery(con, query), error = function(e) return(TRUE))
	if(is(t, "logical")){
		err <- dbGetException(con)
		message(paste("Inserting", disease, "into database caused the following error:", err$errorNum, err$errorMsg, sep=" "))
		return(FALSE)
	} else {
		message(paste("Inserted", disease, "successfully into database", sep=" "))
		return(TRUE)
	}
}

mapDisease <- function(con=NULL, disease){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	disease.db <- getDisease(con)
	map <- match(disease, disease.db[,"name"])
	if(any(is.na(map))){
		new.disease <- unique(disease[which(is.na(map))])
		new.disease <- ifelse(length(new.disease) > 1, paste(new.disease, collapse="','"), new.disease)
		message(paste("Inserting '", new.disease, "' into the DISEASE Table", sep=""))
		if(insertDisease(con, new.disease)){
			disease.db <- getDisease(con)
			map <- match(disease, disease.db[, "name"])
		} else {
			stop(paste("There was an error inserting '", new.disease, "' into the database", sep=" "))
		}
	}
	disease.id <- disease.db[, "id"][map]
	return(disease.id)
}

getTissue <- function(con=NULL){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	tissue <- dbGetQuery(con, "SELECT * FROM TISSUE")
	return(tissue)
}

insertTissue <- function(con=NULL, tissue){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	query <- paste("INSERT INTO TISSUE (name) VALUES ('", tissue, "')", sep="")
	t <- tryCatch(dbSendQuery(con, query), error = function(e) return(TRUE))
	if(is(t, "logical")){
		err <- dbGetException(con)
		message(paste("Inserting", tissue, "into database caused the following error:", err$errorNum, err$errorMsg, sep=" "))
		return(FALSE)
	} else {
		message(paste("Inserted", tissue, "successfully into database", sep=" "))
		return(TRUE)
	}
}

mapTissue <- function(con=NULL, tissue){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	tissue.db <- getTissue(con)
	map <- match(tissue, tissue.db[,"name"])
	if(any(is.na(map))){
		new.tissue <- unique(tissue[which(is.na(map))])
		new.tissue <- ifelse(length(new.tissue) > 1, paste(new.tissue, collapse="','"), new.tissue)
		message(paste("Inserting '", new.tissue, "' into the TISSUE Table", sep=""))
		if(insertTissue(con, new.tissue)){
			tissue.db <- getTissue(con)
			map <- match(tissue, tissue.db[, "name"])
		} else {
			stop(paste("There was an error inserting '", new.tissue, "' into the database", sep=" "))
		}
	}
	tissue.id <- tissue.db[, "id"][map]
	return(tissue.id)
}

getHistology <- function(con=NULL, tissue, disease){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	query <- paste("SELECT * FROM HISTOLOGY WHERE tissue IN ('", paste(tissue, collapse="','"), "') AND disease IN ('", paste(disease, collapse="','"), "')", sep="")
	histology <- dbGetQuery(con, query)
	return(histology)
}

insertHistology <- function(con=NULL, histology){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	query <- paste("INSERT INTO HISTOLOGY (name, tissue, disease) VALUES ", histology, sep="")
	t <- tryCatch(dbSendQuery(con, query), error = function(e) return(TRUE))
	if(is(t, "logical")){
		err <- dbGetException(con)
		message(paste("Inserting", histology, "into database caused the following error:", err$errorNum, err$errorMsg, sep=" "))
		return(FALSE)
	} else {
		message(paste("Inserted", histology, "successfully into database", sep=" "))
		return(TRUE)
	}
}

mapHistology <- function(con=NULL, mappings){
	require("RMySQL")
	if(is.null(con)) stop("Please provide a Database connection object: see ?dbConnect")
	stopifnot(c("tissue", "diseaseabr", "histology") %in% colnames(mappings))
	tissue <- mapTissue(con, mappings$tissue)
	disease <- mapDisease(con, mappings$diseaseabr)
	histology <- data.frame("histology" = mappings$histology, tissue, disease, stringsAsFactors=F)
	histology.db <- getHistology(con, tissue, disease)
	histology.str <- apply(histology, 1, function(x){paste("('", paste(x, collapse="','"), "')", sep="")})
	histology.db.str <- apply(histology.db[ , c("name", "tissue", "disease")], 1, function(x){paste("('", paste(x, collapse="','"), "')", sep="")})
	map <- match(histology.str, histology.db.str)
	if(any(is.na(map))){
		new.histology <- unique(histology.str[which(is.na(map))])
		new.histology <- ifelse(length(new.histology) > 1, paste(new.histology, collapse=","), new.histology)
		message(paste("Inserting", new.histology, "into the HISTOLOGY Table", sep=" "))
		if(insertHistology(con, new.histology)){
			histology.db <- getHistology(con)
			map <- match(histology, histology.db[, "name"])
		} else {
			stop(paste("There was an error inserting", new.histology, "into the database", sep=" "))
		}
	}
	histology.id <- histology.db[, "id"][map]
	return(histology.id)
}

