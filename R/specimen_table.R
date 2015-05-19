
### ---- specimen-table ----
get_specimen_table <- function(impDB, impAlg) {
    spcmTable <- impDB[nzchar(impDB$Extract), c("UFID", "consensusESU", "Country", "Extract", "uuid_idig")]
    spcmTable <- spcmTable[regmatches(spcmTable$Extract, regexpr("^[^,]+", spcmTable$Extract)) %in% dimnames(impAlg)[[1]], ]
    spcmTable <- spcmTable[-match("S0213", spcmTable$Extract), ]
    stopifnot(all(spcmTable$consensusESU %in% load_esuList()))
    spcmTable <- spcmTable[order(spcmTable$consensusESU), ]
    spcmTable$Country <- iconv(spcmTable$Country, "latin1", "ASCII", "")
    spcmTable$Country <- gsub("Runion", "R\\'{e}union", spcmTable$Country, fixed=TRUE)
    spcmTable$Country <- gsub("Nosy B", "Nosy B\\'{e}", spcmTable$Country, fixed=TRUE)
    spcmTable$UFID <- gsub("_", " ", spcmTable$UFID, fixed=TRUE)
    spcmTable$Extract <- gsub("_", " ", spcmTable$Extract, fixed=TRUE)
    names(spcmTable) <- c("Catalog Number", "ESU (consensus)", "Location", "Specimen Nb.", "iDigBio UUID")

    if (FALSE) {
        ## just in case, keeping the code to generate the LaTeX table
        spcmTable <- xtable(spcmTable,
                            caption=c("Specimen information including location, catalog number, and ESU (consensus)"),
                            label="tab:specimen-table")
        headerTable <- paste("\\hline", paste(names(spcmTable), collapse=" & "), "\\\\ \\hline \\endfirsthead \n",
                             "\\caption{(continued) specimen information} \n")
        print.xtable(spcmTable, tabular.environment="longtable", floating=FALSE,
                     hline.after = c(-1, nrow(spcmTable)),
                     add.to.row = list(pos = list(-1, 0),
                         command = c(headerTable, "\\hline \\endhead \n")),
                     sanitize.text.function=function(x) {x}, caption.placement="top",
                     size = "\\tiny",
                     include.rownames=FALSE)
    }

    spcmTable
}
