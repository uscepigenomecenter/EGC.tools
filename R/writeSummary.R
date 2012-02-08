#write sample-level summary

getSampleSummary <- function(pvals)
{
	sample.summary <- apply(pvals, 2, function(x) {length(which(x > 0.05))})
	return(sample.summary)
}

writeSampleSummary <- function(x, filepath='.'){
	disease <- unique(x$diseaseabr)[1]
	pvals <- pvals(x)
	sample.summary <- getSampleSummary(pvals)
	ss <- paste(paste(names(sample.summary), sample.summary, sep="\t"), collapse="\n")
	header <- paste("TCGA.ID", "No. of Failed Probes (p > 0.05)", sep="\t")
	filename <- paste(disease, "sample", "summary", "txt", sep=".")
	cat(header, ss, file=file.path(filepath, filename), sep="\n")
}

getProbeSummary <- function(pvals)
{
	probe.summary <- apply(pvals, 1, function(x) {length(which(x > 0.05))})
	return(probe.summary)
}

writeProbeSummary <- function(x, filepath='.'){
	disease <- unique(x$diseaseabr)[1]
	pvals <- pvals(x)
	probe.summary <- getProbeSummary(pvals)
	ps <- paste(paste(names(probe.summary), probe.summary, sep="\t"), collapse="\n")
	header <- paste("PROBE.ID", "No. of Failed Samples (p > 0.05)", sep="\t")
	filename <- paste(disease, "probe", "summary", "txt", sep=".")
	cat(header, ps, file=file.path(filepath, filename), sep="\n")
}
