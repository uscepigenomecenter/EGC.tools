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
	ss <- data.frame("TCGA.ID"=x$TCGA.ID, "Barcode"=x$barcode, "Batch"=x$TCGA.BATCH, "Num.Failed.Probes"=sample.summary)
	filename <- paste(disease, "sample", "summary", "txt", sep=".")
	write.table(ss, file=file.path(filepath, filename), row.names=FALSE, quote=FALSE, sep="\t")
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

runSampleQC <- function(x, filepath='.')
{
	disease <- unique(x$diseaseabr)[1]
	pvals <- pvals(x)
	sample.summary <- getSampleSummary(pvals)
	if(any(sample.summary > 10000)){
		sample.summary <- sample.summary[which(sample.summary > 10000)]
		failed <- paste(names(sample.summary), collapse=",")
		index <- which(x$barcode %in% names(sample.summary))
		fs <- data.frame("TCGA.ID"=x$TCGA.ID[index], "Barcode"=x$barcode[index], "Histology"=x$histology[index],
				 "Batch"=x$TCGA.BATCH[index], "Num.Failed.Probes"=sample.summary)
		filename <- paste(disease, "sample", "QC", "txt", sep=".")
		message(paste("Samples", failed, "failed QC and will be removed"))
		message(paste("Writing QC log to", file.path(filepath, filename)))
		write.table(fs, file=file.path(filepath, filename), row.names=FALSE, quote=FALSE, sep="\t")
		return(index)
	} else {
		message("All samples passed QC")
		return(NULL)
	}
}
