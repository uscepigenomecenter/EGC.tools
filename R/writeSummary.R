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

plotSampleSummary <- function(x)
{
	disease <- unique(x$diseaseabr)[1]
	pvals <- pvals(x)
	sample.summary <- getSampleSummary(pvals)
	hist(sample.summary, col="salmon", main="Histogram of No. of Failed Probes per Sample")
}

plotDensities <- function(object, controls=NULL, label='Replicate betas'){ # {{{
  if(!is.null(controls)) object <- object[,controls]
  ds <- apply(betas(object), 2, density, from=0, to=1, na.rm=TRUE)
  ymax <- max(do.call(rbind, lapply(ds, function(x) x[['y']])), na.rm=TRUE)
  plot(ds[[1]], main=label, xlim=c(0, 1), ylim=c(0, ymax), xlab='Methylation')
  for(i in seq_along(ds)) lines(ds[[i]])
} # }}}
