strandMe <- function(x) 
{ # {{{ convert +/- coordinates or F/R indicators to c('+','-','*')
  if(is.numeric(x)) ifelse(sign(x) > 0, '+', ifelse(sign(x) < 0, '-', '*'))
  if(all(x %in% c('F','R'))) ifelse(x == 'F', '+', ifelse(x == 'R', '-', '*'))
} # }}}

df2GR <- function(df, keepColumns=FALSE, ignoreStrand=FALSE) 
{ # {{{ adaptation of a utility function from Kasper Daniel Hansen
  stopifnot(class(df) == "data.frame")
  stopifnot(all(c("start", "end") %in% names(df)))
  stopifnot(any(c("chr", "seqnames") %in% names(df)))
  if("seqnames" %in% names(df)) names(df)[names(df) == "seqnames"] <- "chr"
  if(substr(df$chr, 1, 3)[1] != 'chr') df$chr <- paste('chr', df$chr, sep='')
  if(!ignoreStrand && ("strand" %in% names(df))) {
    if(is.numeric(df$strand)) df$strand <- strandMe(df$strand)
    GR <- with(df, GRanges(chr, IRanges(start=start, end=end), strand=strand))
  } else {
    GR <- with(df, GRanges(chr, IRanges(start=start, end=end)))
  }
  if('name' %in% names(df)) {
    names(GR) <- df$name
    df$name <- NULL
  } else {
    names(GR) <- rownames(df)
  }
  if(keepColumns) {
    skipped = c("rangename","chr","start","end","width","strand")
    elementMetadata(GR) <- as(df[, setdiff(names(df), skipped)], "DataFrame")
  }
  return(GR)
} # }}}
