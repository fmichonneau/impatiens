
## code that was used on 2015-05-13 to add iDigBio's UUID to the list of specimens
## included in the study.
if (FALSE) {
    library(ridigbio)
    library(dplyr)

    impDB <- load_impDB("data/impatiensDB.csv") %>%
      filter(nzchar(Extract))

    isUFID <- grepl("^[0-9]+[a-z]?", impDB$UFID)

    UFID_query <- gsub("(^[0-9]+)[a-z]?", "\\1-echinodermata", impDB$UFID[isUFID])

    UFID_idig <- idig_search(rq = list(recordset = "6bb853ab-e8ea-43b1-bd83-47318fc4c345",
                                 catalognumber = UFID_query))

    src <- data.frame(UFID = impDB$UFID[isUFID], UFID_idig = tolower(UFID_query), stringsAsFactors = FALSE)
    res <- data.frame(uuid_idig = UFID_idig$uuid, UFID_idig = UFID_idig$catalognumber, stringsAsFactors = FALSE)

    tmp_res <- left_join(src, res)

    RES <- left_join(impDB, tmp_res)

    write.csv(RES, file = "data/impatiensDB_withuuid.csv", row.names = FALSE)
}
