

getNbFromFile <- function(file) {
    onlyNb <- function(x) {
        gsub("(.+)\\s(.+)$", "\\2", x)
    }
    logFile <- readLines(con=file)
    seedLine <- grep("random", logFile, value=TRUE, ignore.case=TRUE)
    logLinePS <- grep("path .+ pathLikelihood.delta", logFile, value=TRUE)
    logLineSS <- grep("stepping .+ pathLikelihood.delta", logFile, value=TRUE)
    c("Seed"=as.numeric(onlyNb(seedLine)), "logPS"=as.numeric(onlyNb(logLinePS)),
      "logSS"=as.numeric(onlyNb(logLineSS)))
}

getNb <- function(path) {
    owd <- getwd()
    setwd(path)
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
