##' Returns the taxa from a BEAST input file
##'
##' This function works for BEAST >= 1.7.5
##' @title Get taxa from a BEAST input file
##' @param file an XML file used as input for a BEAST analysis
##' @return a vector of character listing the taxa used in the BEAST analysis.
##' @author Francois Michonneau
##' @export
getTaxaBEASTxml <- function(file) {
    doc <- xmlTreeParse(file, getDTD=FALSE)
    r <- xmlRoot(doc)
    unname(sapply(xmlChildren(r[["taxa"]]), function(x) xmlGetAttr(x, "id")))
}
