
source("~/R-scripts/seqManagement.R")
source("~/R-scripts/fasToPhase.R")
source("~/R-scripts/extToLbl.R")
source("~/R-dev/phylothuria/pkg/R/barMonophyletic.R")

runNJonAlg <- function(folder, pattern, db=impDB) {
### Makes NJ trees saved in individual PDF files from alignments
### it also renames the extractions using the extToLbl function
### folder -- is the folder that contains the alignements, also where
###           the PDFs will be saved
### pattern -- regexp to identify the alignment files, supposed .afa extension
### db -- database that links extraction numbers and specimen info (e.g., impDB)
  owd <- getwd()
  listAlg <- list.files(path=folder, pattern=paste(pattern, "afa$", sep=""))
  setwd(folder)
  pdf(file="tmp.pdf", height=40)
  for (i in 1:length(listAlg)) {
    tmpAlg <- read.dna(file=listAlg[i], format="fasta")
    toKeep <- !apply(sapply(tmpAlg, base.freq), 2, function(x) all(is.na(x)))
    cat(length(toKeep), sum(toKeep), "\n")
    tmpTr <- nj(dist.dna(tmpAlg[toKeep, ]))
    tmpTr$tip.label <- extToLblStr(tmpTr$tip.label, db)
    plot.phylo(ladderize(tmpTr), no.margin=T, cex=.5)
    legend("topright", listAlg[i])
  }
  dev.off()
  setwd(owd)
}

makeTraitFile <- function(alg, output, db=impDB) {
### generate the traits file to be used by *BEAST
### alg -- path to alignment in PHYLIP format
### output -- path and filename of where to write the output
  if (!file.exists(alg)) stop(alg, "doesn't exist")
  algF <- read.dna(file=alg)
  tmpRes <- cbind(traits=db$Extract, species=db$consensusESU)
  tmpRes <- tmpRes[nzchar(tmpRes[,1]), ]
  tmpRes[,1] <- gsub(",.+$", "", tmpRes[,1]) # replace extractions by first one if multiple to follow mergeSeq rule
  tmpRes <- subset(tmpRes, tmpRes[,1] %in% dimnames(algF)[[1]])
  allowedESUs <- c("ESU1", "ESU2", "ESU3", "tiger", "RedSea", "Hawaii", "WA", "EP", "Wpac", "gracilis", "Gala", "tigerRedSea", "Medit")
  notAllowed <- tmpRes[,2][! tmpRes[,2] %in% allowedESUs]
  if (length(notAllowed) > 0) warning(paste(unique(notAllowed), collapse=","), ", not found in DB")
  if (any(is.na(tmpRes[,2]))) warning("some NAs found in DB")
  if (length(!nzchar(tmpRes[,2]))) warning("some ESUs not specified")
  write.table(tmpRes, file=output, row.names=FALSE, quote=FALSE, sep="\t")
  TRUE
}

impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv", stringsAsFactors=FALSE)
impExt <- impDB$Extract[nzchar(impDB$Extract)]
stopifnot(ncol(impDB) > 1)

impESU1 <- subset(impDB, molecularESU == "ESU1" | molecularESU == "Hawaii")
impESU1 <- impESU1$Extract[nzchar(impESU1$Extract)]

### Main analysis -- these commands should produce the final output for the published
###  analyses. If it's not the case, it should go somewhere else.

## Create alignments for individual markers
impAlg <- mergeSeq(impExt, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
                   markers=c("16S", "16Sc", "COI", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a"),
                   gblocks=list("ITS" = "-t=d -b4=5 -b5=a -p=n"), justCheck=F)
concatenateAlignments(pattern="afa$", path="/tmp/seq", output="/tmp/seq/20130923.impatiens.phy",
                      partition="/tmp/seq/20130923.impatiens.part", partition.format="nexus",
                      create.conc=TRUE, colsep="", colw=10000)
alg2nex(file="/tmp/seq/20130923.impatiens.phy", partition.file="/tmp/seq/20130923.impatiens.part")

## unfortunately at this stage, still need manual editing of the partition... compare -orig.nex and .nex

## fasToPhase("~/Documents/Impatiens/000.currentFiles/20130710-224540-c0036.afa") # removed N1167 by hand (too much missing data)
## fasToPhase("~/Documents/Impatiens/000.currentFiles/20130710-224540-c0775.afa")
## fasToPhase("~/Documents/Impatiens/000.currentFiles/20130710-224540-H3a.afa")
## fasToPhase("~/Documents/Impatiens/000.currentFiles/20130710-224540-ITS.afa")
## fasToPhase("~/Documents/Impatiens/000.currentFiles/20130710-224540-LSU.afa")
## setwd("../20130711.phase")

## system("./PHASE -d1 20130710-224540-c0036.inp c0036.out")
## system("./PHASE -d1 20130710-224540-c0775.inp c0775.out")
## system("./PHASE -d1 20130710-224540-H3a.inp H3a.out")
## system("./PHASE -d1 20130710-224540-ITS.inp ITS.out")
## system("./PHASE -d1 20130710-224540-LSU.inp LSU.out")

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
