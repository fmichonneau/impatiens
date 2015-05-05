##' Generates a full label from the extraction number
##'
##' The \sQuote{fields} options defines which attributes of the
##' extractions will be used to form the label. The strings of
##' characters must match the column names of the spreadsheet that
##' contains the information about the species. This function is not
##' intended to be used by itself, but rather to be called by
##' \sQuote{extract_to_label} or more generally to be used in
##' conjunction with \sQuote{sapply}.
##' @title Generate labels from the extraction numbers.
##' @param str The character string to be matched and converted
##' @param db The database in CSV format that link the extraction
##' number (\sQuote{str}) with the sample information.
##' @param fields The sequence of information to be collated to form
##' the label.
##' @return A character string corresponding to the formated label.
##' @author Francois Michonneau
extract_to_string_label <- function(str, db, fields=c("Genus", "Species", "consensusESU",
                                          "Country", "UFID", "Extract")) {
  extNb <- str
  extDb <- db$Extract
  dups <- grep(",", extDb)
  tmpRes <- array(, dim=c(1, ncol(db) + 1))
  if (length(dups) > 0) {
      dupDb <- db[dups, ]
      for (i in 1:length(dups)) {
          allDup <- unlist(strsplit(dupDb$Extract[i], ","))
          matDt <- matrix(rep(dupDb[i, ], length(allDup)), nrow=length(allDup),
                          byrow=TRUE)
          tmpDt <- cbind(allDup, matDt)
          tmpRes <- rbind(tmpRes, tmpDt)
      }
  }
  tmpRes <- tmpRes[-1, ]
  tmpRes[, 12] <- tmpRes[, 1]
  tmpRes <- tmpRes[, -1]
  tmpRes <- data.frame(tmpRes)
  names(tmpRes) <- names(db)
  db <- rbind(db, tmpRes)
  fullLbl <- db[match(extNb, db$Extract), ]
  newLbl <- paste("fullLbl$", fields, sep="")
  newLbl <- sapply(newLbl, function(x) unlist(eval(parse(text=x))))
  newLbl <- paste(newLbl, collapse="_")
  newLbl <- gsub("_{2,}", "_", newLbl)
  newLbl
}

##' Converts tip labels from a tree from the extraction number to a full string.
##'
##' For details see the documentation of extract_to_string_label
##' @title Convert tip labels
##' @param tr The tree for which the tip labels need to be converted.
##' @param db The database that link the information between the
##' extraction numbers and the species they belong to.
##' @param ... Additional arguments to be passed to
##' extract_to_string_label. Particularly useful to control elements being used in
##' the labels.
##' @return a tree in the phylo format with converted labels.
##' @author François Michonneau
extract_to_label <- function(tr, db, ...) {
  extNb <- tr$tip.label
  newLbl <- sapply(extNb, extract_to_string_label, db, ...)
  tr$tip.label <- newLbl
  tr
}
