load_impAlg <- function(filename, ...) {
    impAlg <- ape::read.dna(file = filename, format="sequential", ...)
}
