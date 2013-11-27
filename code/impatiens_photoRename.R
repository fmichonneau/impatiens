

impDB <- read.csv(file="data/impatiensDB.csv", stringsAsFactors=FALSE)

lFiles <- list.files(path="../impatiens_specimenPictures", pattern="jpg$|JPG$")

renameByESU <- function(lFiles, impDB) {

    lF <- sapply(lFiles, function(x) {
        unlist(strsplit(x, "-"))[1]
    })

    ufid <- impDB[match(lF, impDB$UFID), ]$consensusESU
    cbind(lFiles, paste(ufid, lFiles, sep="-"))
}

nFiles <- renameByESU(lFiles, impDB)
apply(nFiles, 1, function(x) file.rename(x[1], x[2]))
