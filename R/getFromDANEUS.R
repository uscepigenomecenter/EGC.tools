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
	map <- match(disease, disease.db[,1])
	if(any(is.na(map))){
		new.disease <- unique(disease[which(is.na(map))])
		new.disease <- ifelse(length(new.disease) > 1, paste(new.disease, collapse="','"), new.disease)
		message(paste("Inserting '", new.disease, "' into the DISEASE Table", sep=""))
		if(insertDisease(con, new.disease)){
			disease.db <- getDisease(con)
			map <- match(disease, disease.db[,1])
		} else {
			stop(paste("There was an error inserting '", new.disease, "' into the database", sep=" "))
		}
	}
	disease.id <- disease.db[, 2][map]
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
	map <- match(tissue, tissue.db[,1])
	if(any(is.na(map))){
		new.tissue <- unique(tissue[which(is.na(map))])
		new.tissue <- ifelse(length(new.tissue) > 1, paste(new.tissue, collapse="','"), new.tissue)
		message(paste("Inserting '", new.tissue, "' into the TISSUE Table", sep=""))
		if(insertTissue(con, new.tissue)){
			tissue.db <- getTissue(con)
			map <- match(tissue, tissue.db[,1])
		} else {
			stop(paste("There was an error inserting '", new.tissue, "' into the database", sep=" "))
		}
	}
	tissue.id <- tissue.db[, 2][map]
	return(tissue.id)
}

