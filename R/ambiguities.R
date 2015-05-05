get_alg_files <- function(pattern) {
    list.files(pattern = pattern, path = "data/alignments",  full.names = TRUE)
}

alg_ambiguities <- function(alg_files) {
    res <- sapply(alg_files,  chopper::checkAmbiguity, quiet = TRUE)
    names(res) <- gsub("[0-9]{8}-[0-9]{6}-(.+)\\.afa$", "\\1", basename(alg_files))
    res
}

nAmbPositions_locus <- function(locus, alg_ambiguities) {
    length(unique(unlist(alg_ambiguities[[locus]])))
}

count_nAmbPositions <- function(loci, alg_ambiguities) {
    sum(sapply(loci, nAmbPositions_locus, alg_ambiguities))
}

count_nAmbTotal <- function(loci, alg_ambiguities) {
    length(do.call("c", lapply(loci, function(x) unlist(alg_ambiguities[[x]]))))
}

count_totalNucl <- function(alg) {
    stopifnot(!inherits(alg, "character"))
    sum(! alg %in% c("-",  "?"))
}
