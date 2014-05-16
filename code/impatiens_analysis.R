
### ---- init-phylo ----
setwd("~/Documents/Impatiens/impatiens_phylogeography/")
impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv", stringsAsFactors=FALSE)
library(seqManagement)
library(wesanderson)
source("~/R-scripts/fasToPhase.R")
source("~/R-scripts/extToLbl.R")
source("~/R-dev/phylothuria/pkg/R/barMonophyletic.R")
impExt <- impDB$Extract[nzchar(impDB$Extract)]
stopifnot(ncol(impDB) > 1)



### full impatiens Tree
### ---- impatiens-tree ----
impTree <- read.tree(file="data/allImpatiens.phy")
impTree <- drop.tip(impTree, "S0213") # remove outgroup
impTree <- drop.tip(impTree, "N0057") # remove this one, as only H3a and messes up monophyly, NEED TO REDO analysis!
impTree <- drop.tip(impTree, "G0108") # remove this one, as only mtData and messes up monophyly, NEED TO REDO analysis!
impTree <- ladderize(impTree)
posTips <- treeDepth(impTree) #max(branching.(impTree))
impTree <- extToLbl(impTree, impDB, c("consensusESU", "Country", "UFID", "Extract"))
esuList <- c("ESU1", "ESU2", "ESU3", "gracilis", "tiger", "tigerRedSea", "Medit", "WA",
             "Gala", "EP", "Hawaii", "Wpac", "RedSea")
##esuCol <- hcl(h=seq(0, 360, length=length(esuList)), c=100, l=50)
esuCol <- c("#FE0835", "#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")
impTree$root.edge <- 8 - treeDepth(impTree)
impNodLbl <- as.numeric(impTree$node.label)
impNodLbl[impNodLbl > 1] <- 0.001 # fix little glitch in format conversion, small BEAST posteriors are expressed in scientific notations and only the part before the E is converted which leads to posterior > 1. All are converted to small values .001
impNodLblCol <- rep(NULL, length(impNodLbl))
impNodLblCol[impNodLbl == 1] <- "black"
impNodLblCol[impNodLbl >= .975 & impNodLbl < 1] <- "red"
impNodLblCol[impNodLbl >= .90 & impNodLbl < .975] <- "orange"
plot.phylo(impTree, root.edge=TRUE, show.tip.label=FALSE, x.lim=c(0,10))
nodelabels(text=rep("", length(impNodLblCol)), frame="circ", col=impNodLblCol, bg=impNodLblCol, fg=impNodLblCol, cex=.3)
barMonophyletic(groupLabel=esuList, groupMatch=paste("^", esuList, "_", sep=""), impTree, cex.text=.4,
                cex.plot=.3, extra.space=.5, text.offset=1.02, seg.col=esuCol)
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
### ---- generate-gmyc-trees ----
gmycFactors <- expand.grid(c("strict", "relaxed"),
                           c("yule", "coalexp", "coalcst"))
gmycFactors <- apply(gmycFactors, 1, paste0, collapse="_")
gmycFactors <- c(paste0("20140422.allSeq_", gmycFactors),
                 paste0("20140514.noDup_", gmycFactors))
treeannotatorCmd <-
    c("~/Software/BEASTv1.8.0/bin/./treeannotator -heights ca -burnin 1000")
pathResults <- "~/Documents/Impatiens"
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
    ## each element of the list returns a vector of length 3: mean, range
    list(singleMeanCI=c(tmpSingle$entity[which.max(tmpSingle$likelihood)], ciSingle),
         multiMeanCI=c(tmpMulti$entity[which.max(tmpMulti$likelihood)], ciMulti))
})
names(gmycSumm) <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(gmycSumm))
gmycSumm <- t(data.frame(gmycSumm))
gmycSumm <- data.frame(gmycSumm)
names(gmycSumm) <- c("mean", "low", "high")
tt <- rownames(gmycSumm)
tmpSM <- sapply(tt, function(x) unlist(strsplit(x, "\\."))[2])
tmpSM <- gsub("MeanCI", "", tmpSM)
tmpFac <- strsplit(gsub("\\..+$", "", tt), "_")
tmpSeq <- sapply(tmpFac, function(x) x[1])
tmpClo <- sapply(tmpFac, function(x) x[2])
tmpDem <- sapply(tmpFac, function(x) x[3])
gmycSumm <- cbind(sequences = tmpSeq, clock = tmpClo, demographic = tmpDem,
            analysisType = tmpSM, gmycSumm)

### ---- gmyc-coi-plot ----
levels(gmycSumm$sequences)[levels(gmycSumm$sequences) == "allSeq"] <- "All haplotypes"
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
## library(pegas)
## impHaplo <- subset(impDB, molecularESU %in% c("ESU1"))$Extract
## impHaplo <- impHaplo[nzchar(impHaplo)]
## impCOI <- mergeSeq(impHaplo, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
##                    markers="COI")
## file.copy("/tmp/seq/20121016-114036-COI.afa", "20121016-COI.afa")
## impCOI <- read.dna(file="20121016-COI.afa", format="fasta")


## h <- haplotype(impCOI[-1, ])
## net <- haploNet(h)
## plot(net, size = attr(net, "freq"), scale.ratio = 4, cex = 0.6, labels=F)


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
