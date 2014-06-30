
### ---- init-phylo ----
setwd("~/Documents/Impatiens/impatiens_phylogeography/")
impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv",
                  stringsAsFactors=FALSE)
library(seqManagement)
library(wesanderson)
library(ggplot2)
library(ape)
source("~/R-scripts/fasToPhase.R")
source("code/extToLbl.R")
source("~/R-dev/phylothuria/pkg/R/barMonophyletic.R")
impExt <- impDB$Extract[nzchar(impDB$Extract)]
stopifnot(ncol(impDB) > 1)

## Other common variables
gmycFactors <- expand.grid(c("strict", "relaxed"),
                           c("yule", "coalexp", "coalcst"))
gmycFactors <- apply(gmycFactors, 1, paste0, collapse="_")
gmycFactors <- c(paste0("20140422.allSeq_", gmycFactors),
                 paste0("20140514.noDup_", gmycFactors))
pathResults <- "~/Documents/Impatiens"


### ---- loci-characteristics-data ----
library(seqManagement)
library(pegas)
library(phyloch)
library(xtable)

locFiles <- cutAlignment(algfile="data/20130923.impatiens.phy",
                         partfile="data/20130923.partition-raxml-perlocus",
                         formatin="sequential", format="sequential")

tmpDir <- tempdir()
dir.create(file.path(tmpDir, "data"), recursive=TRUE)

lociChar <- sapply(locFiles, function(fnm) {
    tmpfnm <- file.path(tmpDir, fnm)
    convertGaps(fnm, output=tmpfnm, formatin="sequential", colw=10000)
    removeEmptySeqs(tmpfnm, formatin="sequential", formatout="sequential", gap="-", overwrite=TRUE, colw=10000)
    algChar <- read.dna(file=tmpfnm, format="sequential", as.character=TRUE)
    alg <- read.dna(file=tmpfnm, format="sequential")
    lSeqNonAlign <- apply(algChar, 1, function(x) {
        seq <- paste(as.character(x), collapse="", sep="")
        length(gregexpr("[actgACTG]", seq)[[1]])
    })
    nSeq <- nrow(alg)               # number of individuals sequenced
    nUniq <- nrow(haplotype(alg))   # number of haplotypes/unique sequences
    lSeq <- ncol(alg)               # length of the aligned sequence
    rglSeq <- range(lSeqNonAlign)   # range of raw sequence length
    nSeg <- length(seg.sites(alg))  # number of segregating/variable sites
    nPis <- pis(alg)                # number of parsimony informative sites
    c(nSeq, nUniq, lSeq, rglSeq[2], nSeg, nPis)
})

lociChar <- as.matrix(lociChar)

locNm <- gsub(".+_(.+)\\.phy", "\\1", colnames(lociChar))
colnames(lociChar) <- locNm
rownames(lociChar) <- c("Nind",
                        "K",
                        "bp_alg",
                        "bp_unalg",
                        "S", "S_i")

### ---- loci-characteristics-table ----
lociCharTable <- lociChar
lociCharTable["bp_alg", ] <- paste(lociCharTable["bp_alg", ], " (", lociCharTable["bp_unalg", ], ")", sep="")
lociCharTable <- lociCharTable[- match("bp_unalg", rownames(lociCharTable)), ]

rownames(lociCharTable) <- c("$N$",
                        "$K$",
                        "bp",
                        "$S$",
                        "$S_{i}$")

lociCharTable <- lociCharTable[, c("16S", "COI", "ATP6", "c0036", "c0775", "H3a", "ITS", "LSU")]

locTable <- print(xtable(lociCharTable,
                         caption=c(paste("Characteristics of the loci used for the",
                             "phylogenetic analyses. $N$: number of individuals sequenced,",
                             "$K$: number of haplotypes, bp: length of the aligned (and",
                             "unaligned) sequences, $S$: number of segregating sites,",
                             "$S_{i}$: number of parsimony informative sites.",
                             "The statistics given for ITS are for the ones used in the",
                             "analysis (i.e., after using Gblock)."), # long
                             "Loci characteristics"),
                         label="tab:loci-characteristics"), # short
                  caption.placement="top",
                  sanitize.text.function = function(x) {x},
                  print.results=FALSE)

multiColStr <- paste("\\1 \\& \\\\multicolumn{3}{c}{mtDNA}",
                    "\\\\multicolumn{5}{c}{nucDNA} \\\\\\\\ \\\n",
                    "\\\\cline{2-4} \\\\cline{5-9}", sep="")

cat(gsub("(}\\\n)(\\s+\\\\\\hline)", multiColStr, locTable))

### ---- loci-coverage-data ----
locCov <- sapply(locFiles, function(fnm) {
    tmpfnm <- file.path(tmpDir, fnm)
    convertGaps(fnm, output=tmpfnm, formatin="sequential", colw=10000)
    removeEmptySeqs(tmpfnm, formatin="sequential", formatout="sequential", gap="-", overwrite=TRUE, colw=10000)
    alg <- read.dna(file=tmpfnm, format="sequential")
    dimnames(alg)[[1]]
})
names(locCov) <- gsub(".+_(.+)\\.phy", "\\1", names(locCov))

locWhich <- mapply(function(x, y) cbind(x, rep(y, length(x))), locCov, names(locCov))
locWhich <- do.call("rbind", locWhich)

locWhich <- data.frame(locWhich, stringsAsFactors=FALSE)
names(locWhich) <- c("Extract", "Locus")

locPos <- data.frame("Locus" = unique(locWhich$Locus)[c(1,5,2,3,4,6,7,8)],
                     "begin" = 1:length(unique(locWhich$Locus)),
                     "end"   = 1:length(unique(locWhich$Locus))+1)

locWhich <- merge(locWhich, locPos)

tmpDB <- impDB[, c("Extract", "consensusESU")]
tmpDB <- tmpDB[nzchar(tmpDB$Extract), ]
tmpDB$Extract <- sapply(tmpDB$Extract, function(x) unlist(strsplit(x, ","))[1])

locWhich <- merge(locWhich, tmpDB)

locRank <- data.frame(table(locWhich$Extract))
names(locRank) <- c("Extract", "nLoci")
locRank$Rank <- rank(locRank$nLoci, ties.method="random")

locWhich <- merge(locWhich, locRank)

hasMt <- locWhich[locWhich$Locus %in% c("16S", "COI", "ATP6"), ]

hasMtNuc <- locWhich[locWhich$Extract %in% hasmt$Extract &
                     locWhich$Locus %in% c("c0036", "c0775", "ITS", "LSU", "H3a"), ]


### ---- loci-coverage-plot ----
ggplot(data=locWhich) + geom_segment(aes(x=begin, xend=end,
                            y=Rank, yend=Rank, colour=consensusESU),
                        lineend="round",
                        size=I(1.5)) +
    scale_x_continuous(breaks=1:length(unique(locWhich$Locus))+0.5, labels=locPos$Locus) +
    scale_y_discrete(labels=element_blank()) + ylab("Individuals") + xlab("Loci") +
    theme(legend.position=c(.75,.22),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_blank()) 

### ---- number-ambiguities ----
library(seqManagement)
ambc0036 <- checkAmbiguity(file="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/20130923-185350-c0036.afa",
                           quiet=T)
ambc0775 <- checkAmbiguity(file="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/20130923-185350-c0775.afa",
                           quiet=T)
ambH3a <- checkAmbiguity(file="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/20130923-185350-H3a.afa",
                           quiet=T)
ambITS <- checkAmbiguity(file="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/20130923-185350-ITS.afa",
                           quiet=T)
ambLSU <- checkAmbiguity(file="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/20130923-185350-LSU.afa",
                           quiet=T)
fullAlg <- read.dna(file="data/20130923.impatiens.phy", as.character=TRUE)
nAmbc0036 <- length(unique(unlist(ambc0036)))
nAmbc0775 <- length(unique(unlist(ambc0775)))
nAmbH3a <- length(unique(unlist(ambH3a)))
nAmbITS <- length(unique(unlist(ambITS)))
nAmbLSU <- length(unique(unlist(ambLSU)))
nAmbPositions <- nAmbc0036 + nAmbc0775 + nAmbH3a + nAmbITS + nAmbLSU
nAmbTotal <- length(c(unlist(ambc0036), unlist(ambc0775), unlist(ambH3a), unlist(ambITS), unlist(ambLSU)))
totalNucl <- sum(fullAlg %in% c("a", "c", "t", "g")) 

### ---- impatiens-tree-stats -----
impTree <- read.nexus(file="../20140519.allImpatiens/allimpatiens_strict.tree.nex")

### full impatiens Tree
### ---- impatiens-tree ----
library(ape)
source("code/getPosteriors.R")
source("code/extToLbl.R")
impTree <- read.beast(file="data/20140519.impTree.nex")
if( !file.exists("data/20140519.impTree-beast.phy") )
    write.tree(impTree, file="data/20140519.impTree-beast.phy")
if (! file.exists("data/RAxML_bootstrap_nooutgroup.phy")) {
    allRaxmlBS <- read.tree(file="data/RAxML_bootstrap.result")
    noOutgroup <- lapply(allRaxmlBS, function(tr) drop.tip(tr, "S0213"))
    class(noOutgroup) <- "multiPhylo"
    write.tree(noOutgroup, file="data/RAxML_bootstrap_nooutgroup.phy")
}
if (length(list.files(pattern="annotateBEASTtree$", path="data/")) != 3) {
    rxmlCmd <- paste("~/Software/RAxML-8.0.1/./raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA",
                     "-p 12345 -f b -t data/20140519.impTree-beast.phy",
                     "-z data/RAxML_bootstrap_nooutgroup.phy -T8 -n annotateBEASTtree")
    system(rxmlCmd)
    system("mv *.annotateBEASTtree data/")
}
## noOutgroup <- read.tree(file="data/RAxML_bootstrap_nooutgroup.phy")
## ppNoOutGroup <- prop.part(noOutgroup)
## pcNoOutGroup <- prop.clades(impTree, part=ppNoOutGroup)
## save(pcNoOutGroup, file="data/pcNoOutGroup.RData")
load("data/pcNoOutGroup.RData")
impTree <- ladderize(impTree)
posTips <- max(branching.times(impTree))
impTree <- extToLbl(impTree, impDB, c("consensusESU", "Country", "UFID", "Extract"))
esuList <- c("ESU1", "ESU2", "ESU3", "gracilis", "tiger", "tigerRedSea", "Medit", "WA",
             "Gala", "EP", "Hawaii", "Wpac", "RedSea")
##esuCol <- hcl(h=seq(0, 360, length=length(esuList)), c=sample(c(50, 100), size=length(esuList), replace=TRUE), l=50)
##esuCol <- c("#FE0835", "#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")
esuCol <- c(
"#25162E",
"#4F261D",
"#7D4954",
"#D6E4DB",
"#8FA6C4",
"#9E8498",
"#3C366A",
"#B7565A",
"#D0AD87",
"#EED780",
"#B7CF64",
"#6F8042",
"#576266")

esuCol <- c("#51574a", "#447c69", "#74c493", "#8e8c6d", "#e4bf80", "#e9d78e", "#e2975d",
            "#f19670", "#e16552", "#c94a53", "#be5168", "#a34974", "#993767", "#65387d",
            "#4e2472", "#9163b6", "#e279a3", "#e0598b", "#7c9fb0", "#5698c4", "#9abf8e")


#impTree$root.edge <- 8 - treeDepth(impTree)

#impNodLbl[impNodLbl > 1] <- 0.001 # fix little glitch in format conversion, small BEAST posteriors are expressed in scientific notations and only the part before the E is converted which leads to posterior > 1. All are converted to small values .001
impNodLblCol <- rep(NULL, length(impNodLbl))
impNodLblCol[impNodLbl == 1] <- "black"
impNodLblCol[impNodLbl >= .99 & impNodLbl < 1] <- "red"
impNodLblCol[impNodLbl >= .90 & impNodLbl < .99] <- "orange"

impNodLblTxt <- rep("", length(impNodLblCol))
impAllNds <- 1:length(impNodLbl)
impKeepNds <- impAllNds[!is.na(impNodLblCol)]

plot.phylo(impTree, root.edge=TRUE, show.tip.label=FALSE)#, x.lim=c(0,10))
nodelabels(text=rep("", length(impKeepNds)), node=impKeepNds+Ntip(impTree),
               frame="circ", col=impNodLblCol[!is.na(impNodLblCol)],
               bg=impNodLblCol[!is.na(impNodLblCol)],
               fg=impNodLblCol[!is.na(impNodLblCol)],
               cex=.3)

barMonophyletic(groupLabel=esuList, groupMatch=paste("^", esuList, "_", sep=""), impTree, cex.text=.4,
                cex.plot=.3, extra.space=.5, text.offset=1.02,
                seg.col=rev(esuCol))
axis(side=1, at=0:8, labels=paste("-", 8:0, sep=""))
legend(x=0, y=50, pch=16, col=c("black", "red", "orange"), legend=c("$PP = 1$", "$0.975 \\leq PP < 1$", "$ 0.9 \\leq PP < 0.975$"))


### StarBeast summary results (marginal likelihood summaries)
### ---- starbeast-summary ----
sbeastOrig <- read.csv(file="data/starbeastResults.csv")
sbeast <- sbeastOrig[, c("groupings", "runs", "chainLength", "dataIncluded",
                         "PS_logLik", "SS_logLik")]
sbeast <- sbeast[grep("[0-9]$", sbeast$runs), ] # runs with X at the ends are not meant to be included
wMtSB <- subset(sbeast, chainLength == "long" & dataIncluded == "all")
bppsMat <- matrix(wMtSB[, c("PS_logLik")], nrow=2); bppsMat <- bppsMat[, -8]
bppsMat <- bppsMat - min(bppsMat[,1]) #bppsMat[, 1]
bpssMat <- matrix(wMtSB[, c("SS_logLik")], nrow=2); bpssMat <- bpssMat[, -8]
bpssMat <- bpssMat - min(bpssMat[,1]) #bpssMat[, 1]
lbls <- unique(wMtSB$groupings); lbls <- lbls[-length(lbls)]
noCOISB <- subset(sbeast, chainLength == "long" & dataIncluded == "noCOI")
noCOIbpssMat <- matrix(noCOISB[, c("SS_logLik")], nrow=2)
noCOIbpssMat <- noCOIbpssMat - max(noCOIbpssMat[, 1])
noCOIlbls <- unique(noCOISB$groupings)
par(mfrow=c(1, 2))
bpss <- barplot(bpssMat, ylab="Difference in log-marginal likelihoods",
                beside=T, main="With COI", ylim=c(-50, 5))
mtext(side=1, at=colMeans(bpss), line=0, text=as.character(lbls), las=2, adj=0)
text(bpss[1:2], bpssMat[,1] + c(-1.5, 1.5), c("*", "*"))
abline(h=-5, lty=2)
noCOIbpss <- barplot(noCOIbpssMat, beside=T, main="Without COI",
                     ylim=c(-50, 5))
mtext(side=1, at=colMeans(noCOIbpss), line=0, text=as.character(noCOIlbls),
      las=2, adj=0)
abline(h=-5, lty=2)
text(noCOIbpss[1:2], noCOIbpssMat[,1] + c(-1.5, 1.5), c("*", "*"))
## noMtbppsMat <- matrix(noMtSBeast[, c("PS_logLik")], nrow=2); 
## noMtbppsMat <- noMtbppsMat - min(noMtbppsMat[, 1])
## noMtbpssMat <- matrix(noMtSBeast[, c("SS_logLik")], nrow=2); 
## noMtbpssMat <- noMtbpssMat - min(noMtbpssMat[, 1])
## noMtlbls <- unique(noMtSBeast$groupings); 
## noMtbpss <- barplot(noMtbpssMat, ylab="Difference in log-marginal likelihoods",
##                     beside=T, main="Stepping-stone sampling", ylim=c(-5,5))
## mtext(side=1, at=colMeans(noMtbpss), line=0, text=noMtlbls, las=2, adj=0)


### ---- sm-starbeast-summary ----
bpps <- barplot(bppsMat, ylab="Difference in log-marginal likelihoods",
                beside=T, #main="Path sampling",
                ylim=c(-50, 0))
mtext(side=1, at=colMeans(bpps), line=0, text=as.character(lbls), las=2, adj=0)
text(bpps[1:2], bppsMat[,1] - 1.5, c("*", "*"))
abline(h=-5, lty=2)
## noMtbpps <- barplot(noMtbppsMat, ylab="Difference in log-marginal likelihoods",
##                     beside=T, main="Path sampling", ylim=c(-5,5))
## mtext(side=1, at=colMeans(noMtbpps), line=0, text=noMtlbls, las=2, adj=0)

### NJ analyses
### ---- nj-coi ----
source("~/R-scripts/findGroups.R")
impCOI <- read.nexus.data(file="data/20131009_impatiens_COIuniq.nex")
impNJ <- nj(dist.dna(as.DNAbin(impCOI)))
impNJ$edge.length[impNJ$edge.length < 0] <- 0
impNJ <- root(impNJ, "S0213", resolve.root=T)
impNJ <- drop.tip(impNJ, "S0213")
trNJ4 <- as(ladderize(impNJ), "phylo4")
grpImp <- findGroups(trNJ4, threshold=.015)
tipLabels(grpImp) <- paste(tipData(grpImp)$Group, tipLabels(grpImp), sep="_")
grpLbl <- paste("^", 1:max(tipData(grpImp)), "_", sep="")
par(mai=c(0.5,0,2,0), xpd=T)
plot(as(grpImp, "phylo"), cex=.5, show.tip.label=T, no.margin=F, label.offset=0)
barMonophyletic(1:max(tipData(grpImp)), grpLbl, as(grpImp, "phylo"), cex.plot=.5,
                cex.text=.5, extra.space=0.01, bar.at.tips=TRUE,
                include.tip.label=TRUE)
add.scale.bar()

### GMYC analyses
### ---- gmyc-target-tree    ----
library(XML)
xmlBeastGetTaxa <- function(file) {
    doc <- xmlTreeParse(file, getDTD=FALSE)
    r <- xmlRoot(doc)
    unname(sapply(xmlChildren(r[["taxa"]]), function(x) xmlGetAttr(x, "id")))
}
targetTreeFull <- read.nexus(file="~/Documents/Impatiens/20140519.allImpatiens/allimpatiens_strict.tree.nex")
labelsCoi <- xmlBeastGetTaxa(file="~/Documents/Impatiens/20140422.allSeq_relaxed_yule/20140422.allSeq_relaxed_yule.xml")
targetTreeLblAll <- targetTreeFull$tip.label
targetTree <- drop.tip(targetTreeFull, targetTreeLblAll[! targetTreeLblAll %in% labelsCoi])
write.tree(targetTree, file="data/20140616.impatiens_targetTree.tre")

### ---- gmyc-allmt-analyses ----
##  dir.create("~/Documents/Impatiens/20140627.impatiens_allMt")
## concatenateAlignments(pattern="20130923-.*(16Sc?|COI|ATP6)\\.afa$",
##                       path="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/",
##                       output="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.phy",
##                       partition="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.part",
##                       partition.format="nexus",
##                       format="seq", colw=1000, colsep="")
## mtImp <- read.dna(file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.phy",
##                   format="seq")
mtImp <- mtImp[-match("S0213", dimnames(mtImp)[[1]]), ]
write.dna(mtImp, file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt_noS0213.phy",
          format="seq", colw=10000, colsep="")
alg2nex(file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt_noS0213.phy",
        partition.file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.part")

### ---- generate-gmyc-trees ----
## pathResults, gmycFactors are defined in init-phylo
treeannotatorCmd <-
    c("~/Software/BEASTv1.8.0/bin/./treeannotator -heights ca -burnin 1000")
inFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".trees"))
outFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".tre.nex"))
## can't parallize this, too much RAM needed.
## foreach (i = 1:length(gmycFactors)) %dopar% {
for (i in 1:length(gmycFactors)) {
    if (file.exists(inFile[i])) {
        tmpCmd <- paste(treeannotatorCmd, inFile[i], outFile[i])
        if (! file.exists(outFile[i])) {
            system(tmpCmd)
        }
        else {
            message(outFile[i], " already exists.")
        }
    }
    else {
        message(inFile[i], " doesn't exist.")
    }
}

treeannotatorCmdWithTarget <- paste(
    c("~/Software/BEASTv1.8.0/bin/./treeannotator -heights ca -burnin 1000",
      "-target ~/Documents/Impatiens/impatiens_phylogeography/data/20140616.impatiens_targetTree.tre"),
    collapse=" ")
inFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".trees"))
outFileTarget <- file.path(pathResults, gmycFactors, paste0(gmycFactors, "_withTarget.tre.nex"))
## can't parallize this, too much RAM needed.
## foreach (i = 1:length(gmycFactors)) %dopar% {
for (i in 1:length(gmycFactors)) {
    if (file.exists(inFile[i])) {
        tmpCmd <- paste(treeannotatorCmdWithTarget, inFile[i], outFileTarget[i])
        if (! file.exists(outFileTarget[i])) {
            system(tmpCmd)
        }
        else {
            message(outFileTarget[i], " already exists.")
        }
    }
    else {
        message(inFileTarget[i], " doesn't exist.")
    }
}

### ---- gmyc-coi-results ----
### depends on previous chunk to generate the trees, but once they are there
###   no need to be run again
## library(splits)
## library(doMC)
## library(ggplot2)
## registerDoMC(cores=7)
## gmycRes <- foreach (i = 1:length(outFile)) %dopar% {
##     if (! file.exists(outFile[i])) {
##         message(outFile[i], " doesn't exist.")
##     }
##     else {
##         tmpTr <- read.nexus(file=outFile[i])
##         list(simpleGmyc = gmyc(tmpTr), multiGmyc = gmyc(tmpTr, method="m"))
##     }
## }
## names(gmycRes) <- outFile
## save(gmycRes, file="data/gmycRes.RData")
load("data/gmycRes.RData")

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

 ## ggplot(data=gmycSumm, aes(x=ageTree, y=mean, colour=interaction(sequences, clock, demographic),
 ##           shape=analysisType)) + geom_point()

### ---- gmyc-coi-plot ----
levels(gmycSumm$sequences)[levels(gmycSumm$sequences) == "allSeq"] <- "All Sequences"
levels(gmycSumm$sequences)[levels(gmycSumm$sequences) == "noDup"]  <- "Unique haplotypes"
levels(gmycSumm$clock)[levels(gmycSumm$clock) == "strict"] <- "Strict clock"
levels(gmycSumm$clock)[levels(gmycSumm$clock) == "relaxed"] <- "Relaxed clock"
ggplot(data=gmycSumm, aes(x=demographic, y=mean,
           ymin=low, ymax=high, color=analysisType)) +
    geom_linerange(width=.2, linetype=2, position=position_dodge(width = 0.6)) +
    geom_point(position=position_dodge(width = 0.6)) + 
    ylab("Estimated number of species") +
    scale_y_continuous(breaks=seq(12,30,by=2)) +
    facet_grid( ~ sequences + clock) +
    labs(x = "Tree priors") + 
    theme(legend.position="top", legend.title=element_blank(),
          panel.background = element_rect(fill = "gray95"),
          panel.grid.major = element_line(colour = "white", size=0.1),
          panel.grid.minor = element_line(NA)) +
    scale_color_manual(values = wes.palette(5, "Zissou")[c(1, 5)],
                        labels=c("multi-threshold GMYC", "single threshold GMYC"))

### ---- gmyc-tree-plot ----
## find the tree with the most conservative estimate
nspecies <- sapply(gmycRes, function(x) {
    tmpSingle <- x$simpleGmyc
    tmpSingle$entity[which.max(tmpSingle$likelihood)]
})
resKeep <- gmycRes[[which.min(whichKeep)]]$simpleGmyc
stopifnot(identical(class(resKeep), "gmyc"))
keepTree <- resKeep$tree
## Get threshold
keepThreshold <- resKeep$threshold.time[which.max(resKeep$likelihood)]
## Find edge length for edge leading to node associate with threshold
##   so we can draw line in the middle of the edge
keepEdgeWhich <- names(branching.times(keepTree)[which(branching.times(keepTree) == -keepThreshold)])
keepEdgeLength <- keepTree$edge.length[which(keepTree$edge[, 2] == as.integer(keepEdgeWhich))]
plot(ladderize(keepTree), cex=.5)
abline(v=max(branching.times(keepTree))+keepThreshold-(keepEdgeLength/2), col="red", lwd=2)

### ---- median-tmrca-WAgroup ----
## Get mean and median for node corresponding to MRCA for all
##  individuals found on each side of isthmus of Panama as estimated
##  with various parameters used to test GMYC.  Prior on this node is
##  log normal prior with a log(mean) of 1.5, a log(standard
##  deviation) of 0.75 and an offset of 2.5. This distribution
##  translates into a median age of divergence set at 6.98 millions
##  years (MY) ago, (95% confidence of [3.53, 21.99])

## logFiles <- file.path(pathResults, gmycFactors, paste(gmycFactors, ".log", sep=""))
## stopifnot(all(file.exists(logFiles)))
## resTmrca <- lapply(logFiles, function(logF) {
##     tmpDt <- read.table(logF, header=TRUE, colClasses="numeric")
##     stopifnot(nrow(tmpDt) == 10001)
##     colid <- grep("tmrca", names(tmpDt))
##     stopifnot(length(colid) == 1)
##     rg <- 1001:nrow(tmpDt)
##     list(mean=mean(tmpDt[rg, colid]),
##          median=median(tmpDt[rg, colid]))
## })
## names(resTmrca) <- logFiles
## save(resTmrca, file="data/resTmrca.RData")
load(file="data/resTmrca.RData")
meanTmrca <- sapply(resTmrca, function(x) x$mean)
medianTmrca <- sapply(resTmrca, function(x) x$median)

### starBEAST PPS
### ---- starbeast-pps ----
library(starbeastPPS)

impxml <- read.starbeast(beast.xml="/home/francois/Photos/impatiens_analyses/000.allESU1/starbeastPPS_allESU1_longSampling_run1/20131230.allESU1_longSampling.xml", "/home/francois/Photos/impatiens_analyses/000.allESU1/starbeastPPS_allESU1_longSampling_run1/combined")

impcoal <- analyze.coalescent(impxml, "~/Software/msdir/")

impseq <- analyze.sequences(impxml, "/home/francois/Software/Seq-Gen.v1.3.3")

### bGMYC analysis
### ---- bgmyc ----
library(bGMYC)

resSingle <- bgmyc.singlephy(trBeast, mcmc=50000, burnin=1, thinning=10, t1=2, t2=100, start=c(1,1,25))

### ---- ClockStaR ----
## system("cd ~/R-dev/; git clone git@github.com:sebastianduchene/ClockstaR.git")
## install("~/R-dev/ClockstaR")
library(ClockstaR2)
optim.trees.interactive()

#####################

### Just for ESU1
## impESU1Alg <- mergeSeq(impESU1, output="/tmp", seqFolder="~/Documents/seqRepository",
##                        markers=c("16Sc", "16S", "COI", "ATP6", "H3a", "LSU"))
## concatenateAlignments(patter="afa$", path="/tmp", output="/tmp/20130202.impESU1.phy",
##                       partition="/tmp/20130202.impESU1.part", partition.format="nexus",
##                       create.conc=T, colsep="", colw=10000)


### to get all extractions that have at least another marker other than COI
## impAlg <- mergeSeq(impExt, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
##                    markers=c("16S", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a", "COI"),
##                    justCheck=T)
## impAlg <- impAlg[(impAlg$"16S" | impAlg$ATP6 | impAlg$c0036 | impAlg$ITS | impAlg$LSU | impAlg$H3a), ]
## impExt2 <- rownames(impAlg)
## impAlg <-  mergeSeq(impExt2, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
##                    markers=c("16S", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a", "COI"))
## concatenateAlignments(patter="afa$", path="/tmp/seq", output="/tmp/20130411.impatiens.phy",
##                       partition="/tmp/20130411.impatiens.part", partition.format="nexus",
##                       create.conc=T, colsep="", colw=10000)

## ### compare individual topologies (NJ based)
## trc0036 <- read.tree(file="20130507-c0036.tre")
## trc0775 <- read.tree(file="20130507-c0775.tre")
## trLSU <- read.tree(file="20130507-LSU.tre")
## trITS <- read.tree(file="20130507-ITS.tre")
## trH3a <- read.tree(file="../20130507-H3a.tre")

## trc0036 <- extToLbl(trc0036, impDB)
## trc0775 <- extToLbl(trc0775, impDB)
## trLSU <- extToLbl(trLSU, impDB)
## trITS <- extToLbl(trITS, impDB)
## trH3a <- extToLbl(trH3a, impDB)

## pdf(file="/tmp/trees.pdf")
## plot(trc0036, no.margin=T, cex=.5)
## plot(trc0775, no.margin=T, cex=.5)
## plot(trLSU, no.margin=T, cex=.5)
## plot(trITS, no.margin=T, cex=.5)
## plot(H3a, no.margin=T, cex=.5)
## dev.off()


### run RAxML to get reduced file
## phy2nex(file="~/Documents/Impatiens/20130507.impatiens_partfinder/20130507.impatiens.phy",
##         partition.file="~/Documents/Impatiens/20130507.impatiens_partfinder/20130507.impatiens.part")

## rownames(impAlg) <- extToLblStr(rownames(impAlg), impDB)

## tmpTr <- read.tree(file="/tmp/seq/RAxML_bestTree.T1")
## tmpTr <- root(tmpTr, grep("S0213", tmpTr$tip.label))
## tmpTr$tip.label <- extToLblStr(tmpTr$tip.label, impDB)
## tmpTr <- ladderize(tmpTr)
## posMarkers <- seq(1, by=.02, length.out=ncol(impAlg))
## tipOrder <- tmpTr$tip.label[tmpTr$edge[which(tmpTr$edge %in% 1:Ntip(tmpTr))]]
## pdf(file="/tmp/20121022.tree.pdf", height=20)
## plot(tmpTr, no.margin=TRUE, cex=.4)
## points(rep(posMarkers[1], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder, ]$"16Sc", 15, 22), cex=.6)
## points(rep(posMarkers[2], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"16S", 15, 22), cex=.6)
## points(rep(posMarkers[3], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"COI", 15, 22), cex=.6)
## points(rep(posMarkers[4], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"ATP6", 15, 22), cex=.6)
## points(rep(posMarkers[5], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"H3a", 15, 22), cex=.6)
## points(rep(posMarkers[6], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"LSU", 15, 22), cex=.6)
## text(posMarkers, rep(Ntip(tmpTr)+1, length(posMarkers)), c("16Sc", "16S", "COI", "ATP6", "H3a", "LSU"), srt=90, cex=.4)
## dev.off()


### summary table
## impAlgNm <- extToLblStr(rownames(impAlg), impDB)
## impAlg$ESU <- sapply(impAlgNm, function(x) unlist(strsplit(x, "_"))[3])

## library(reshape)
## impSum <- melt(impAlg, id.vars="ESU", variable_name="marker")

## xtabs(~ ESU + marker, data=impSum, subset=value)

### haplotype networks
 library(pegas)
esu1HaploData <- subset(impDB, consensusESU %in% c("ESU1"))$Extract
esu1HaploData <- esu1HaploData[nzchar(esu1HaploData)]
esu1HaploData <- gsub(",.+$", "", esu1HaploData)

esu2HaploData <- subset(impDB, consensusESU %in% c("ESU2"))$Extract
esu2HaploData <- esu2HaploData[nzchar(esu2HaploData)]
esu2HaploData <- gsub(",.+$", "", esu2HaploData)

impCOI <- read.dna(file="data/20130923.impatiens_COI.phy", format="seq")

esu1Seq <- impCOI[match(esu1HaploData, dimnames(impCOI)[[1]]), ]
esu1ToRm <- which(sapply(esu1Seq, function(x) any(is.na(base.freq(x)))))
esu1Seq <- esu1Seq[-esu1ToRm, ]

esu2Seq <- impCOI[match(esu2HaploData, dimnames(impCOI)[[1]]), ]
esu2ToRm <- which(sapply(esu2Seq, function(x) any(is.na(base.freq(x)))))
esu2Seq <- esu2Seq[-esu2ToRm, ]


esu1Haplotypes <- haplotype(esu1Seq)
esu1HaploNet <- haploNet(esu1Haplotypes)
plot(esu1HaploNet, size=attr(esu1HaploNet, "freq"), bg=ziss, scale.ratio=20, cex=.6)


esu2Haplotypes <- haplotype(esu2Seq)
esu2HaploNet <- haploNet(esu2Haplotypes)
plot(esu2HaploNet, size=attr(esu2HaploNet, "freq"))

### test GMYC
## library(splits)
## trBeast <- read.nexus(file="/home/francois/Documents/Impatiens/20130201.beast_COI_PB_strictClock/20130201-COI-PB-strict.tre")
## singleGmycBeast <- gmyc(trBeast)
## multiGmycBeast <- gmyc(trBeast, method="multiple")

## trNJ <- read.tree(file="~/Documents/Impatiens/20121105.impatiens_beast/20130131.COI-NJ_ultra.nex")
## trNJ <- drop.tip(trNJ, "S0213")
## trNJ <- multi2di(trNJ)
## singleGmycNJ <- gmyc(trNJ)
## multiGmycNJ <- gmyc(trNJ, method="multiple")

## trBeastCOI <- read.nexus(file="~/Documents/Impatiens/20130131.impatiens_beast_COI/20130131-111957-COI_tree.nex")
## trBeastCOI <- drop.tip(trBeastCOI, "S0213")
## singleGmycBeastCOI <- gmyc(trBeastCOI)
## multiGmycBeastCOI <- gmyc(trBeastCOI, method="multiple")

#####################
