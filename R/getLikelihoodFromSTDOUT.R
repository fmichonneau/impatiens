

getNbFromFile <- function(file, fast=TRUE) {
    onlyNb <- function(x) {
        gsub("(.+)\\s(.+)$", "\\2", x)
    }
    if (fast) { ## this option doesn't return the seed
        logFile <- system(paste("tail -300", file), intern=TRUE)
    }
    else {
        logFile <- readLines(con=file)
    }
    seedLine <- grep("random", logFile, value=TRUE, ignore.case=TRUE)
    logLinePS <- grep("path .+ pathLikelihood.delta", logFile, value=TRUE)
    logLineSS <- grep("stepping .+ pathLikelihood.delta", logFile, value=TRUE)
    c("Seed"=as.numeric(onlyNb(seedLine)),
      "logPS"=as.numeric(onlyNb(logLinePS)),
      "logSS"=as.numeric(onlyNb(logLineSS)))
}

getNb <- function(path) {
    owd <- setwd(path)
    allStdOut <- list.files(path=path, recursive=TRUE, pattern="STDOUT")
    res <- vector("list", length(allStdOut))
    for (i in 1:length(allStdOut)) {
        res[[i]] <- getNbFromFile(allStdOut[i])
    }
    setwd(owd)
    fNm <- gsub("(.+)/(.+)/(.+)$", "\\1", allStdOut)
    names(res) <- fNm
    res
}
