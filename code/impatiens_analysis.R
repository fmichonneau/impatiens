
### ---- init-phylo ----
setwd("~/Documents/Impatiens/impatiens_phylogeography/")
impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv",
                  stringsAsFactors=FALSE)
library(seqManagement)
library(ggplot2)
library(ape)
library(phylobase)

## external scripts (need to be moved)
source("~/R-dev/phylothuria/pkg/R/barMonophyletic.R")
source("~/R-scripts/fasToPhase.R")

## scripts
source("code/gmycGroups.R")
source("code/multiplot.R")
source("code/extToLbl.R")
source("code/getPosteriors.R")

## quick check on data
impExt <- impDB$Extract[nzchar(impDB$Extract)]
stopifnot(ncol(impDB) > 1)

## beast tree
## TODO -- change path
impTree <- read.beast(file="../20140519.allImpatiens/allimpatiens_strict.tree.nex")
impTree4 <- as(impTree, "phylo4")

## DNA alignment for all data
impAlg <- ape::read.dna(file="data/20130923.impatiens.phy", format="seq")

## ESU list
esuList <- c("ESU1", "ESU2", "ESU3", "gracilis", "tiger", "tigerRedSea", "Medit", "WA",
             "Gala", "EP", "Hawaii", "Wpac", "RedSea")

## Other common variables
gmycFactors <- expand.grid(c("strict", "relaxed"),
                           c("yule", "coalexp", "coalcst"))
gmycFactors <- apply(gmycFactors, 1, paste0, collapse="_")
gmycFactors <- c(paste0("20140422.allSeq_", gmycFactors),
                 paste0("20140514.noDup_", gmycFactors))
pathResults <- "~/Documents/Impatiens"

## colors
impPal <- c(
    "#89B151", # tenenbaum bright green = Medit
    "#AF420A", # life aquatic dark orange = ESU3_Deep
    "#FD6467", # grand budapest hotel pink = ESU3_PNG
    "#01abe9", # life aquatic bright blue = ESU1
    "#1b346c", # life aquatic dark blue = ESU3
    "#f54b1a", # life aquatic bright orange = gracilis
    "#e5c39e", # life aquatic tan = Hawaii
    "#c3ced0", # life aquatic light blue = RedSea
    "#EBCC2A", # life aquatic yellow = Wpac
    "#446455", # tenenbaum dark green = Gala
    "#4CBDD5", # darjeeling light blue = WA
    "#FF0000", # darjeeling red = tigerRedSea
    "#F8B0BB", # pink = tiger
    "#B6F9B1", # another green = ESU1_Lizard
    "#FFF196", # mr fox yellow = EP
    "#8061AD" # rushmore purple = ESU2
    #"#020202"  # not quite black
    )
names(impPal) <- c("Medit", "ESU3_Deep", "ESU3_PNG", "ESU1", "ESU3", "gracilis",
                   "Hawaii", "RedSea", "Wpac", "Gala", "WA", "tigerRedSea", "tiger",
                   "ESU1_Lizard", "EP", "ESU2")

## image(1:length(impPal), 1, as.matrix(1:length(impPal)), col=impPal, xlab="impPal",
##        ylab="", xaxt="n", yaxt="n", bty="n")


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
                             "$K$: number of unique sequences, bp: length of the aligned",
                             "(unaligned) sequences, $S$: number of segregating sites,",
                             "$S_{i}$: number of parsimony informative sites.",
                             "The statistics given for ITS are for the ones used in the",
                             "analysis (i.e., after using Gblock)."), # long
                             "Loci characteristics"),  # short
                         label="tab:loci-characteristics"),
                  caption.placement="top",
                  sanitize.text.function = function(x) {x},
                  print.results=FALSE)

multiColStr <- paste("\\1 \\& \\\\multicolumn{3}{c}{mtDNA} \\& ",
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
    scale_colour_manual(values=impPal) +
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


### full impatiens Tree
### ---- impatiens-tree ----
library(ape)
source("code/extToLbl.R")
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

impNodLbl <- impTree$posterior
impNodLblCol <- rep(NULL, length(impNodLbl))
impNodLblCol[impNodLbl >= .99] <- "black"

impAllNds <- 1:length(impNodLbl)
impKeepNds <- impAllNds[!is.na(impNodLblCol)]

par(mai=c(1,0,0,0))
plot.phylo(impTree, root.edge=TRUE, show.tip.label=FALSE, x.lim=c(0,12))
nodelabels(text=rep("", length(impKeepNds)), node=impKeepNds+Ntip(impTree),
               frame="circ", col=impNodLblCol[!is.na(impNodLblCol)],
               bg=impNodLblCol[!is.na(impNodLblCol)],
               ##fg=impNodLblCol[!is.na(impNodLblCol)],
               cex=.3)
barMonophyletic(groupLabel=esuList, groupMatch=paste("^", esuList, "_", sep=""), impTree, cex.text=.4,
                cex.plot=.3, extra.space=.1, text.offset=1.02,
                seg.col=impPal[esuList])
abline(v=max(branching.times(impTree)) - c(1:4, 6:10), lty=1, col="gray70")
abline(v=max(branching.times(impTree) - 5), lty=2, col="gray50")
axis(side=1, at=max(branching.times(impTree)) - 0:11, labels=c(0, paste("-", 1:10, sep=""), NA))
legend(x=0, y=10, pch=16, col="black", legend=c("PP $\\geq$ 0.99"), bg="white")

### RAxML tree
### ---- raxml-tree ----
raxmlImpTr <- read.tree(file="data/raxml_fullcomplex/RAxML_bipartitions.result",)
raxmlImpTr <- drop.tip(root(raxmlImpTr, "S0213"), "S0213")
raxmlImpTr <- ladderize(raxmlImpTr)

raxmlImpTr <- extToLbl(raxmlImpTr, impDB, c("consensusESU", "Country", "UFID", "Extract"))

ndCol <- character(length(raxmlImpTr$node.label))
ndCol[as.numeric(raxmlImpTr$node.label) >= 90] <- "black"
ndCol[as.numeric(raxmlImpTr$node.label) >= 80 &
      as.numeric(raxmlImpTr$node.label) < 90] <- "red"
ndAll <- 1:length(ndCol)
ndKeep <- ndAll[nzchar(ndCol)]

par(mai=c(0,0,0,0))
plot(raxmlImpTr, cex=.6, show.tip.label=FALSE, x.lim=c(0,.245))
barMonophyletic(groupLabel=esuList, groupMatch=paste("^", esuList, "_", sep=""),
                raxmlImpTr, cex.text=.4, cex.plot=.3, extra.space=.155,
                text.offset=1.02, seg.col=impPal[esuList])
nodelabels(text=rep("", length(ndKeep)), node=ndKeep+Ntip(raxmlImpTr),
           bg=ndCol[nzchar(ndCol)], frame="circ", cex=.275)
add.scale.bar()


### Per locus trees
### ---- per-locus-trees ----
loci <- c("c0036", "c0775", "H3a", "mtDNA", "rDNA")
raxFiles <- file.path("data/raxml_perlocus", paste0("RAxML_bipartitions.", loci))
stopifnot(all(sapply(raxFiles, file.exists)))
uniqExt <- regmatches(impDB$Extract, regexpr("^[^,]+", impDB$Extract))

par(mai=c(0,0,.5,0))
layout(matrix(c(4,1,2,4,3,5), 2, 3, byrow=TRUE))
for (i in 1:length(raxFiles)) {
    raxmlTr <- read.tree(file=raxFiles[i])
    raxmlCol <- impDB[match(raxmlTr$tip.label, uniqExt), "consensusESU"]
    raxmlCol <- impPal[raxmlCol]    
    plotTr <- ladderize(raxmlTr)
    plotTr$tip.label <- gsub("_", " ", plotTr$tip.label)    
    ndCol <- character(length(plotTr$node.label))
    ndCol[as.numeric(plotTr$node.label) >= 80] <- "black"
    ndAll <- 1:length(ndCol)
    ndKeep <- ndAll[nzchar(ndCol)]
    plot(plotTr, tip.color=raxmlCol, main=loci[i], cex=.5)
    nodelabels(text=rep("", length(ndKeep)), node=ndKeep+Ntip(plotTr), bg=ndCol[nzchar(ndCol)], frame="circ",
               cex=.3)
    add.scale.bar()
    if (i == 4) {
        textLgd <- c("Medit", "WA", "Gala", "EP", "ESU2", "tiger", "tigerRedSea", "ESU3", "gracilis",
                     "RedSea", "ESU1", "Wpac", "Hawaii")
        legend(x=0, y=70, legend=textLgd, col=impPal[textLgd], lty=1, lwd=3)
    }
}

### StarBeast summary results (marginal likelihood summaries)
### ---- starbeast-summary ----
sbeastOrig <- read.csv(file="data/starbeastResults.csv")
sbeast <- sbeastOrig[, c("groupings", "runs", "chainLength", "dataIncluded",
                         "PS_logLik", "SS_logLik")]
sbeast <- sbeast[grep("[0-9]$", sbeast$runs), ] # runs with X at the ends are not meant to be included

## wMtSB <- subset(sbeast, chainLength == "long" & dataIncluded == "all")
## bppsMat <- matrix(wMtSB[, c("PS_logLik")], nrow=2); bppsMat <- bppsMat[, -ncol(bppsMat)]
## bppsMat <- bppsMat - min(bppsMat[,1]) #bppsMat[, 1]
## bpssMat <- matrix(wMtSB[, c("SS_logLik")], nrow=2); bpssMat <- bpssMat[, -ncol(bpssMat)]
## bpssMat <- bpssMat - min(bpssMat[,1]) #bpssMat[, 1]
## lbls <- unique(wMtSB$groupings); lbls <- lbls[-length(lbls)]
## noCOISB <- subset(sbeast, chainLength == "long" & dataIncluded == "noCOI")
## noCOIbpssMat <- matrix(noCOISB[, c("SS_logLik")], nrow=2)
## noCOIbpssMat <- noCOIbpssMat - max(noCOIbpssMat[, 1])
## noCOIlbls <- unique(noCOISB$groupings)

## par(mfrow=c(1, 2))
## bpss <- barplot(bpssMat, ylab="Difference in log-marginal likelihoods",
##                 beside=T, main="With COI", ylim=c(-50, 5))
## mtext(side=1, at=colMeans(bpss), line=0, text=as.character(lbls), las=2, adj=0)
## text(bpss[1:2], bpssMat[,1] + c(-1.5, 1.5), c("*", "*"))
## abline(h=-5, lty=2)
## noCOIbpss <- barplot(noCOIbpssMat, beside=T, main="Without COI",
##                      ylim=c(-50, 5))
## mtext(side=1, at=colMeans(noCOIbpss), line=0, text=as.character(noCOIlbls),
##       las=2, adj=0)
## abline(h=-5, lty=2)
## text(noCOIbpss[1:2], noCOIbpssMat[,1] + c(-1.5, 1.5), c("*", "*"))

sbCol <- wes.palette(5, "Zissou")[c(1, 5)]

sbSummAll <- subset(sbeast, chainLength == "long" & dataIncluded == "all" & groupings != "random")
meanESU1SS <- mean(subset(sbSummAll, groupings == "allESU1")$SS_logLik)
meanESU1PS <- mean(subset(sbSummAll, groupings == "allESU1")$PS_logLik) 
sbSummAll$stdSS <- sbSummAll$SS_logLik - meanESU1SS
sbSummAll$stdPS <- sbSummAll$PS_logLik - meanESU1PS
sbSummAll$dataIncluded <- "All data"
sbSummNoCOI <- subset(sbeast, chainLength == "long" & dataIncluded == "noCOI" & groupings != "random")
meanESU1SSnoCOI <- mean(subset(sbSummNoCOI, groupings == "allESU1")$SS_logLik)
meanESU1PSnoCOI <- mean(subset(sbSummNoCOI, groupings == "allESU1")$PS_logLik) 
sbSummNoCOI$stdSS <- sbSummNoCOI$SS_logLik - meanESU1SSnoCOI
sbSummNoCOI$stdPS <- sbSummNoCOI$PS_logLik - meanESU1SSnoCOI
sbSummNoCOI$dataIncluded <- "No COI"

BFHawaii <- round(2 * mean(subset(sbSummAll, groupings == "noHawaii")$stdPS), 1)
BFWpac <- round(2 * mean(subset(sbSummAll, groupings == "noWpac")$stdPS), 1)
BFRedSea <- round(2 * mean(subset(sbSummAll, groupings == "noRedSea")$stdPS), 1)
BFsplit  <- round(2 * mean(subset(sbSummAll, groupings == "allESU1split")$stdPS), 1)

ggplot(sbSummAll, aes(x=groupings, y=stdSS)) + geom_point(position="dodge", colour=sbCol[1]) +
    stat_summary(fun.y = mean, geom="point", colour=sbCol[2], size=3) +
    geom_hline(yintercept=-5, width=.2, col=sbCol[2], linetype=2) +
    theme(legend.position="top", legend.title=element_blank(),
          panel.background = element_rect(fill = "gray95"),
          panel.grid.major = element_line(colour = "white", size=0.1),
          panel.grid.minor = element_line(NA)) +
    labs(y="Difference in Marginal log-Likelihood", x="Models") +
    facet_grid( ~ dataIncluded)


## sbNoCoiPlotSS <- ggplot(sbSummNoCOI, aes(x=groupings, y=stdSS)) + geom_point(position="dodge", colour=sbCol[1]) +
##     stat_summary(fun.y = mean, geom="point", colour=sbCol[2], size=3) +
##     geom_hline(yintercept=-5, width=.2, col=sbCol[2], linetype=2) +
##     theme(legend.position="top", legend.title=element_blank(),
##           panel.background = element_rect(fill = "gray95"),
##           panel.grid.major = element_line(colour = "white", size=0.1),
##           panel.grid.minor = element_line(NA)) +
##     labs(y="", x="Models") +
##     facet_grid(~ dataIncluded)


#multiplot(sbAllPlotSS, sbNoCoiPlotSS, layout=matrix(c(1,1,2), nrow=1))


### ---- sm-starbeast-summary ----
## bpps <- barplot(bppsMat, ylab="Difference in log-marginal likelihoods",
##                 beside=T, #main="Path sampling",
##                 ylim=c(-50, 0))
## mtext(side=1, at=colMeans(bpps), line=0, text=as.character(lbls), las=2, adj=0)
## text(bpps[1:2], bppsMat[,1] - 1.5, c("*", "*"))
## abline(h=-5, lty=2)


ggplot(sbSummAll, aes(x=groupings, y=stdPS)) + #, ymin=0, ymax=stdPS)) +
    geom_point(aes(x=groupings, y=stdPS), colour=sbCol[1], size=1) +
    stat_summary(fun.y = mean, geom="point", colour=sbCol[2], size=3) +
    geom_hline(yintercept=-5, width=.2, col=sbCol[2], linetype=2) +
    theme(legend.position="top", legend.title=element_blank(),
          panel.background = element_rect(fill = "gray95"),
          panel.grid.major = element_line(colour = "white", size=0.1),
          panel.grid.minor = element_line(NA)) +
    labs(y="Difference in Marginal log-Likelihood", x="Models") +
    facet_grid( ~ dataIncluded)        

## sbNoCoiPlotPS <- ggplot(sbSummNoCOI, aes(x=groupings, y=stdPS)) + geom_point(position="dodge", colour=sbCol[1]) +
##     stat_summary(fun.y = mean, geom="point", colour=sbCol[2], size=3) +
##     geom_hline(yintercept=-5, width=.2, col=sbCol[2], linetype=2) +
##     theme(legend.position="top", legend.title=element_blank(),
##           panel.background = element_rect(fill = "gray95"),
##           panel.grid.major = element_line(colour = "white", size=0.1),
##           panel.grid.minor = element_line(NA)) +
##     labs(y="", x="Models") +
##     facet_grid(~ dataIncluded)

## multiplot(sbAllPlotSS, sbNoCoiPlotSS, layout=matrix(c(1,1,2), nrow=1))
                    

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

### ---- gmyc-tree-plot-prep ----
## find the tree with the most conservative estimate
load("data/gmycRes.RData")
nspecies <- sapply(gmycRes, function(x) {
    tmpSingle <- x$simpleGmyc
    tmpSingle$entity[which.max(tmpSingle$likelihood)]
})
nmTrees <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(gmycRes))
nmTrees <- gsub("^allSeq_", "All sequences, ", nmTrees)
nmTrees <- gsub("^noDup_", "Haplotypes, ", nmTrees)
nmTrees <- gsub("relaxed_", "Relaxed Clock, ", nmTrees)
nmTrees <- gsub("strict_", "Strict Clock, ", nmTrees)
nmTrees <- gsub("yule$", "Yule", nmTrees)
nmTrees <- gsub("coalexp", "Coalescent (Exponential growth)", nmTrees)
nmTrees <- gsub("coalcst", "Coalescent (Constant)", nmTrees)

allSeqGmEd <- lapply(gmycRes[1:6], function(x) gmycEdges(x$simpleGmyc))
uniqSeqGmEd <- lapply(gmycRes[7:12], function(x) gmycEdges(x$simpleGmyc))

### ---- gmyc-tree-plot-allSeq ---- 
par(mfrow=c(3, 2))
for (i in 1:6) {
    resKeep <- gmycRes[[i]]$simpleGmyc
    stopifnot(identical(class(resKeep), "gmyc"))
    keepTree <- resKeep$tree
    ## Get threshold
    keepThreshold <- resKeep$threshold.time[which.max(resKeep$likelihood)]
    ## Find edge length for edge leading to node associate with threshold
    ##   so we can draw line in the middle of the edge
    keepBrTimes <- sort(branching.times(keepTree))
    keepWhichEdge <- which(keepBrTimes == -keepThreshold)
    ## We place the threshold line between where the thershold point is and
    ##  the following (in branching order) node in the tree
    thresLine <- max(branching.times(keepTree)) - 
        mean(c(keepBrTimes[keepWhichEdge], keepBrTimes[keepWhichEdge + 1]))
    plotTree <- ladderize(keepTree)
    colEdges <- rep("black", nrow(plotTree$edge))
    edgesToChange <- gmycMatchEdges(allSeqGmEd[[i]], plotTree)
    nSpp <- length(edgesToChange)
    if (nSpp == 14) {
        colToUse <- impPal[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3", "gracilis",
                             "tiger", "Gala", "Medit", "ESU3_Deep", "ESU1_Lizard", "RedSea",
                             "Wpac")]
    } else if (nSpp == 16) {
        colToUse <- impPal[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3",
                             "gracilis", "tiger", "Gala", "Medit", "ESU3_PNG",
                             "ESU1_Lizard", "RedSea", "Wpac", "ESU3_Deep",
                             "tigerRedSea")]
    }
    else stop("Houston, we have a problem.")
    for (j in 1:nSpp) {
        colEdges[edgesToChange[[j]]] <- colToUse[j]
    }
    par(mai= c(1, 0, 0, 0))
    plot(plotTree, show.tip.label=FALSE, no.margin=TRUE,
         x.lim=c( 2 * thresLine - max(keepBrTimes), max(keepBrTimes)),
         edge.color=colEdges, edge.width=0.5)
    mtext(paste(LETTERS[i], nmTrees[i], sep=". "), side=1, line=-1, cex=.7)
    segments(x0=thresLine, x1=thresLine, y0=10,
             y1=Ntip(keepTree), col="red", lwd=2, lty=2)
}

### ---- gmyc-tree-plot-uniqSeq ---- 
par(mfrow=c(3, 2))
for (i in 7:12) {
    k <- i - 6
    resKeep <- gmycRes[[i]]$simpleGmyc
    stopifnot(identical(class(resKeep), "gmyc"))
    keepTree <- resKeep$tree
    ## Get threshold
    keepThreshold <- resKeep$threshold.time[which.max(resKeep$likelihood)]
    ## Find edge length for edge leading to node associate with threshold
    ##   so we can draw line in the middle of the edge
    keepBrTimes <- sort(branching.times(keepTree))
    keepWhichEdge <- which(keepBrTimes == -keepThreshold)
    ## We place the threshold line between where the thershold point is and
    ##  the following (in branching order) node in the tree
    thresLine <- max(branching.times(keepTree)) - 
        mean(c(keepBrTimes[keepWhichEdge], keepBrTimes[keepWhichEdge + 1]))
    plotTree <- ladderize(keepTree)
    colEdges <- rep("black", nrow(plotTree$edge))
    edgesToChange <- gmycMatchEdges(uniqSeqGmEd[[k]], plotTree)
    nSpp <- length(edgesToChange)
    if (nSpp == 14) {
        colToUse <- impPal[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3", "gracilis",
                             "tiger", "Gala", "Medit", "ESU3_Deep", "ESU1_Lizard", "RedSea",
                             "Wpac")]
    }  else if (nSpp == 16) {
        colToUse <- impPal[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3",
                             "gracilis", "tiger", "Gala", "Medit", "ESU3_PNG",
                             "ESU1_Lizard", "RedSea", "Wpac", "ESU3_Deep",
                             "tigerRedSea")]
    } else stop("Houston, we have a problem.")
    for (j in 1:nSpp) {
        colEdges[edgesToChange[[j]]] <- colToUse[j]
    }
    par(mai= c(1, 0, 0, 0))
    plot(plotTree, show.tip.label=FALSE, no.margin=TRUE,
         x.lim=c( 2 * thresLine - max(keepBrTimes), max(keepBrTimes)),
         edge.color=colEdges, edge.width=0.5)
    mtext(paste(LETTERS[i], nmTrees[i], sep=". "), side=1, line=-1, cex=.7)
    segments(x0=thresLine, x1=thresLine, y0=10,
             y1=Ntip(keepTree), col="red", lwd=2, lty=2)
}

### ---- specimen-table ----
library(xtable)
spcmTable <- impDB[nzchar(impDB$Extract), c("UFID", "consensusESU", "Country", "Extract")]
spcmTable <- spcmTable[regmatches(spcmTable$Extract, regexpr("^[^,]+", spcmTable$Extract)) %in% dimnames(impAlg)[[1]], ]
spcmTable <- spcmTable[-match("S0213", spcmTable$Extract), ]
stopifnot(all(spcmTable$consensusESU %in% esuList))
spcmTable <- spcmTable[order(spcmTable$consensusESU), ]
spcmTable$Country <- iconv(spcmTable$Country, "latin1", "ASCII", "")
spcmTable$Country <- gsub("Runion", "R\\'{e}union", spcmTable$Country, fixed=TRUE)
spcmTable$Country <- gsub("Nosy B", "Nosy B\\'{e}", spcmTable$Country, fixed=TRUE)
spcmTable$UFID <- gsub("_", " ", spcmTable$UFID, fixed=TRUE)
spcmTable$Extract <- gsub("_", " ", spcmTable$Extract, fixed=TRUE)
names(spcmTable) <- c("Catalog Nb.", "ESU (consensus)", "Location", "Specimen Nb.")
spcmTable <- xtable(spcmTable,
                    caption=c("Specimen information including location, catalog number, and ESU (consensus)",
                        "Specimen information"),
                    label="tab:specimen-table")
headerTable <- paste("\\hline", paste(names(spcmTable), collapse=" & "), "\\\\ \\hline \\endfirsthead \n",
                     "\\caption{(continued) specimen information} \n")
print.xtable(spcmTable, tabular.environment="longtable", floating=FALSE,
             hline.after = c(-1, nrow(spcmTable)),
             add.to.row = list(pos = list(-1, 0),
                 command = c(headerTable, "\\hline \\endhead \n")),
             sanitize.text.function=function(x) {x}, caption.placement="top",
             include.rownames=FALSE)

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
