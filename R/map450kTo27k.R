slots.27k.to.450k <- function(slots) { # {{{
  rs = paste('R0', 1:6, sep='')
  slotses = c(paste(rs, 'C01', sep=''), paste(rs, 'C02', sep=''))
  names(slotses) = toupper(letters[1:length(slotses)])
  return(slotses[toupper(slots)])
} # }}}
