gmycGroups <- function(gmycOut) {

    stopifnot(gmycOut$method == "single")

    tr <- gmycOut$tree
    threshold <- -gmycOut$threshold.time[which.max(gmycOut$likelihood)]
    bt <- branching.times(tr)
    tr4 <- as(tr, "phylo4")

    grps <- integer(length(tr$tip.label))
    names(grps) <- tr$tip.label
    grpNb <- 1L
    treeDepth <- max(branching.times(tr))

    for (i in 1:length(tr$tip.label)) {
        if (grps[i] == 0) {
            tmpAnc <- phylobase::ancestors(tr4, i)
            tmpBt <- sapply(tmpAnc, function(x) nodeDepth(tr4, x))
            tmpBt <- treeDepth - tmpBt
            dThreshold <- tmpBt - threshold
            if (any(dThreshold < 0)) {
                tmpGrp <- phylobase::descendants(tr4,
                                                 tmpAnc[which(dThreshold == max(dThreshold[dThreshold <= 0]))], type="tips")
            }
            else {
                tmpGrp <- getNode(tr4, tipLabels(tr4)[i])
            }
            grps[names(tmpGrp)] <- grpNb
            grpNb <- grpNb + 1L
            if (length(grps) > length(tr$tip.label)) browser()
        }
        else {
            next
        }
    }
    grps
}

gmycEdges <- function(gmycOut) {
    ## returns a list of edges (in phylobase format) associated for which
    ## each element is a groups identified with gmyc
    gmGrp <- gmycGroups(gmycOut)
    tr <- gmycOut$tree
    tr4 <- as(tr, "phylo4")
    nGrp <- max(gmGrp)
    edList <- lapply(1:16, function(x) {
        tips <- names(gmGrp)[gmGrp == x]
        if (length(tips) > 1) {
            anc <- MRCA(tr4, tips)
            desc <- phylobase::descendants(tr4, anc, "all")
            getEdge(tr4, desc)
        }
        else
            getEdge(tr4, tips)
    })
}

gmycMatchEdges <- function(listEdges, tr) {
    ## returns index of each of the edges
    ## listEdges is output of gmycEdges
    trEdges <- paste(tr$edge[,1], tr$edge[,2], sep="-")
    lapply(listEdges, function(e) match(e, trEdges))
}


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
