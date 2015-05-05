load_impDB <- function(filename) {
    impDB <- read.csv(file = filename, stringsAsFactors=FALSE)
    stopifnot(ncol(impDB) > 1) # to make sure the data is read in correctly
    impDB$SingleExtract[nzchar(impDB$Extract)] <- regmatches(impDB$Extract, regexpr("^[^,]+", impDB$Extract))
    impDB
}
