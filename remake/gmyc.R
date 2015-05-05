get_gmyc_summary <- function(gmycRes) {
    gmycSumm <- lapply(gmycRes, function(x) {
                           tmpSingle <- x$simpleGmyc
                           tmpMulti <- x$multiGmyc
                           ciSingle <- range(tmpSingle$entity[tmpSingle$likelihood > (max(tmpSingle$likelihood) - 2)])
                           ciMulti <- range(tmpMulti$entity[tmpMulti$likelihood > (max(tmpMulti$likelihood) - 2)])
                           timeSingle <- tmpSingle$threshold.time[which.max(tmpSingle$likelihood)]
                           ## each element of the list returns a vector of length 4: mean, range, threshold time (not implemented
                           ##   for the multi-threshold)
                           list(singleMeanCITime=c(tmpSingle$entity[which.max(tmpSingle$likelihood)], ciSingle, timeSingle),
                                multiMeanCITime=c(tmpMulti$entity[which.max(tmpMulti$likelihood)], ciMulti, NA))
                       })
    names(gmycSumm) <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(gmycSumm))
    gmycSumm <- t(data.frame(gmycSumm))
    gmycSumm <- data.frame(gmycSumm)
    names(gmycSumm) <- c("mean", "low", "high", "threshold.time")
    tt <- rownames(gmycSumm)
    tmpSM <- sapply(tt, function(x) unlist(strsplit(x, "\\."))[2])
    tmpSM <- gsub("MeanCITime", "", tmpSM)
    tmpFac <- strsplit(gsub("\\..+$", "", tt), "_")
    tmpSeq <- sapply(tmpFac, function(x) x[1])
    tmpClo <- sapply(tmpFac, function(x) x[2])
    tmpDem <- sapply(tmpFac, function(x) x[3])
    gmycTrees <- lapply(gmycRes, function(x) x$simpleGmyc$tree)
    ageTrees <- lapply(gmycTrees, function(x) rep(max(branching.times(x)), 2))
    timeSingle <- lapply(gmycRes, function(x) x$simpleGmyc$threshold.time[which.max(x$simpleGmyc$likelihood)])
    names(ageTrees) <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(ageTrees))
    gmycSumm <- cbind(sequences = tmpSeq, clock = tmpClo, demographic = tmpDem,
                      analysisType = tmpSM, gmycSumm, ageTree = unlist(ageTrees))
    gmycSumm
}
