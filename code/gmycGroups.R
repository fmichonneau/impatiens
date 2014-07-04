
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



if (FALSE) {
    
sGmycRes <- lapply(gmycRes, function(x) x$simpleGmyc)
gmycGrps <- lapply(sGmycRes, function(x) gmycGroups(x))
names(gmycGrps) <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(gmycGrps))
gmycGrpsDt <- data.frame(gmycGrps[1:6])

assgmt <- apply(gmycGrpsDt, 1, function(x) sum(duplicated(x)) != ncol(gmycGrpsDt) - 1)

apply(gmycGrpsDt, 1, function(x) x) #identical(x[1], x))

ordr <- sapply(1:14, function(x) rownames(gmycGrpsDt)[gmycGrpsDt[, 2] == x])
gmycGrpsDt <- gmycGrpsDt[unlist(ordr), ]
gmycGrpsDt$samples <- rownames(gmycGrpsDt)
gmycGrpsDt$yymin <- 1:length(unique(gmycGrpsDt$samples))
##gmycGrpsDt$yymax <- gmycGrpsDt$yymin + 0.5
gmycGrpsDt <- melt(gmycGrpsDt, id.vars=c("samples", "yymin"))
gmycGrpsDt$xxmin <- as.numeric(gmycGrpsDt$variable)
gmycGrpsDt$xxmax <- as.numeric(gmycGrpsDt$variable) + 0.5

ggplot(gmycGrpsDt,
       aes(x=xxmin, xend=xxmax, y=yymin, yend=yymin, colour=factor(value), size=1)) +
    geom_segment() + scale_colour_discrete()

}
