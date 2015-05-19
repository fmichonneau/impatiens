### ---- common-variables ----
## Loci
nuclear_loci <- get_nuclear_loci()
mt_loci <- get_mt_loci()
all_loci <- get_all_loci()

## ESUs
esuList <- load_esuList()

## Colors
impPal <- load_impPal()

## GMYC factors
gmycFactors <- load_gmycFactors()

## Trees
impTree4 <- load_impTree4(impTree)

### ---- ambiguities ----
nAmbPositions <- count_nAmbPositions(loci = nuclear_loci, alg_ambiguities = alg_amb)
nAmbTotal <- count_nAmbTotal(loci = nuclear_loci, alg_ambiguities = alg_amb)
totalNucl <- count_totalNucl(impAlg_char)

### ---- star-beast-stats ----
sBeastESU <- merge(data.frame(SingleExtract=sBeastTaxa), impDB)
sBeastESU <- table(sBeastESU$consensusESU)

### ---- loci-coverage-data ----
hasMt <- locus_graph_data[locus_graph_data$Locus %in% c("16S", "COI", "ATP6"), ]

hasMtNuc <- locus_graph_data[locus_graph_data$Extract %in% hasMt$Extract &
                     locus_graph_data$Locus %in% c("c0036", "c0775", "ITS", "LSU", "H3a"), ]


## ### GMYC analyses
## TODO --- will need to transition all of this to remake
## ### ---- gmyc-target-tree    ----
## targetTreeFull <- read.nexus(file="~/Documents/Impatiens/20140519.allImpatiens/allimpatiens_strict.tree.nex")
## labelsCoi <- xmlBeastGetTaxa(file="~/Documents/Impatiens/20140422.allSeq_relaxed_yule/20140422.allSeq_relaxed_yule.xml")
## targetTreeLblAll <- targetTreeFull$tip.label
## targetTree <- drop.tip(targetTreeFull, targetTreeLblAll[! targetTreeLblAll %in% labelsCoi])
## write.tree(targetTree, file="data/20140616.impatiens_targetTree.tre")

## ### ---- gmyc-allmt-analyses ----
## ##  dir.create("~/Documents/Impatiens/20140627.impatiens_allMt")
## ## concatenateAlignments(pattern="20130923-.*(16Sc?|COI|ATP6)\\.afa$",
## ##                       path="~/Documents/Impatiens/000.impatiens_datafiles/a000.currentFiles/",
## ##                       output="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.phy",
## ##                       partition="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.part",
## ##                       partition.format="nexus",
## ##                       format="seq", colw=1000, colsep="")
## ## mtImp <- read.dna(file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.phy",
## ##                   format="seq")
## mtImp <- mtImp[-match("S0213", dimnames(mtImp)[[1]]), ]
## write.dna(mtImp, file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt_noS0213.phy",
##           format="seq", colw=10000, colsep="")
## alg2nex(file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt_noS0213.phy",
##         partition.file="~/Documents/Impatiens/20140627.impatiens_allMt/20140627.impatiens_allMt.part")

## ### ---- generate-gmyc-trees ----
## ## pathResults, gmycFactors are defined in init-phylo
## treeannotatorCmd <-
##     c("~/Software/BEASTv1.8.0/bin/./treeannotator -heights ca -burnin 1000")
## inFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".trees"))
## outFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".tre.nex"))
## ## can't parallize this, too much RAM needed.
## ## foreach (i = 1:length(gmycFactors)) %dopar% {
## for (i in 1:length(gmycFactors)) {
##     if (file.exists(inFile[i])) {
##         tmpCmd <- paste(treeannotatorCmd, inFile[i], outFile[i])
##         if (! file.exists(outFile[i])) {
##             system(tmpCmd)
##         }
##         else {
##             message(outFile[i], " already exists.")
##         }
##     }
##     else {
##         message(inFile[i], " doesn't exist.")
##     }
## }

## treeannotatorCmdWithTarget <- paste(
##     c("~/Software/BEASTv1.8.0/bin/./treeannotator -heights ca -burnin 1000",
##       "-target ~/Documents/Impatiens/impatiens_phylogeography/data/20140616.impatiens_targetTree.tre"),
##     collapse=" ")
## inFile <- file.path(pathResults, gmycFactors, paste0(gmycFactors, ".trees"))
## outFileTarget <- file.path(pathResults, gmycFactors, paste0(gmycFactors, "_withTarget.tre.nex"))
## ## can't parallize this, too much RAM needed.
## ## foreach (i = 1:length(gmycFactors)) %dopar% {
## for (i in 1:length(gmycFactors)) {
##     if (file.exists(inFile[i])) {
##         tmpCmd <- paste(treeannotatorCmdWithTarget, inFile[i], outFileTarget[i])
##         if (! file.exists(outFileTarget[i])) {
##             system(tmpCmd)
##         }
##         else {
##             message(outFileTarget[i], " already exists.")
##         }
##     }
##     else {
##         message(inFileTarget[i], " doesn't exist.")
##     }
## }

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
## load("data/gmycRes.RData")
## later I converted this to RDS --> data/gmyc_results.rds


## ggplot(data=gmycSumm, aes(x=ageTree, y=mean, colour=interaction(sequences, clock, demographic),
 ##           shape=analysisType)) + geom_point()


### ---- starbeast-summary ----
meanESU1SS <- mean(subset(star_beast_summary_all, groupings == "allESU1")$SS_logLik)

BFHawaii <- round(2 * mean(subset(star_beast_summary_all, groupings == "noHawaii")$stdSS), 1)
BFWpac <- round(2 * mean(subset(star_beast_summary_all, groupings == "noWpac")$stdSS), 1)
BFRedSea <- round(2 * mean(subset(star_beast_summary_all, groupings == "noRedSea")$stdSS), 1)
BFsplit  <- round(2 * mean(subset(star_beast_summary_all, groupings == "allESU1split")$stdSS), 1)
BFrandom <- round(2 * (mean(star_beast_results[star_beast_results$groupings=="random" &
                                                 !is.na(star_beast_results$SS_logLik), "SS_logLik"]) - meanESU1SS), 1)

BFHawaiiNoCOI <- round(2 * mean(subset(star_beast_summary_noCOI, groupings == "noHawaii")$stdSS), 1)
BFWpacNoCOI <- round(2 * mean(subset(star_beast_summary_noCOI, groupings == "noWpac")$stdSS), 1)
BFRedSeaNoCOI <- round(2 * mean(subset(star_beast_summary_noCOI, groupings == "noRedSea")$stdSS), 1)

BFHawaiiNoMt <- round(2 * mean(subset(star_beast_summary_noMt, groupings == "noHawaii")$stdSS), 1)
BFWpacNoMt <- round(2 * mean(subset(star_beast_summary_noMt, groupings == "noWpac")$stdSS), 1)
BFRedSeaNoMt <- round(2 * mean(subset(star_beast_summary_noMt, groupings == "noRedSea")$stdSS), 1)




## ## sbNoCoiPlotSS <- ggplot(star_beast_summary_noCOI, aes(x=groupings, y=stdSS)) + geom_point(position="dodge", colour=sbCol[1]) +
## ##     stat_summary(fun.y = mean, geom="point", colour=sbCol[2], size=3) +
## ##     geom_hline(yintercept=-5, width=.2, col=sbCol[2], linetype=2) +
## ##     theme(legend.position="top", legend.title=element_blank(),
## ##           panel.background = element_rect(fill = "gray95"),
## ##           panel.grid.major = element_line(colour = "white", size=0.1),
## ##           panel.grid.minor = element_line(NA)) +
## ##     labs(y="", x="Models") +
## ##     facet_grid(~ dataIncluded)


## #multiplot(sbAllPlotSS, sbNoCoiPlotSS, layout=matrix(c(1,1,2), nrow=1))


## ### NJ analyses
## ### ---- nj-coi ----
## source("~/R-scripts/findGroups.R")
## impCOI <- read.nexus.data(file="data/20131009_impatiens_COIuniq.nex")
## impNJ <- nj(dist.dna(as.DNAbin(impCOI)))
## impNJ$edge.length[impNJ$edge.length < 0] <- 0
## impNJ <- root(impNJ, "S0213", resolve.root=T)
## impNJ <- drop.tip(impNJ, "S0213")
## trNJ4 <- as(ladderize(impNJ), "phylo4")
## grpImp <- findGroups(trNJ4, threshold=.015)
## tipLabels(grpImp) <- paste(tipData(grpImp)$Group, tipLabels(grpImp), sep="_")
## grpLbl <- paste("^", 1:max(tipData(grpImp)), "_", sep="")
## par(mai=c(0.5,0,2,0), xpd=T)
## plot(as(grpImp, "phylo"), cex=.5, show.tip.label=T, no.margin=F, label.offset=0)
## barMonophyletic(1:max(tipData(grpImp)), grpLbl, as(grpImp, "phylo"), cex.plot=.5,
##                 cex.text=.5, extra.space=0.01, bar.at.tips=TRUE,
##                 include.tip.label=TRUE)
## add.scale.bar()




## ### ---- median-tmrca-WAgroup ----
## ## Get mean and median for node corresponding to MRCA for all
## ##  individuals found on each side of isthmus of Panama as estimated
## ##  with various parameters used to test GMYC.  Prior on this node is
## ##  log normal prior with a log(mean) of 1.5, a log(standard
## ##  deviation) of 0.75 and an offset of 2.5. This distribution
## ##  translates into a median age of divergence set at 6.98 millions
## ##  years (MY) ago, (95% confidence of [3.53, 21.99])

## ## logFiles <- file.path(pathResults, gmycFactors, paste(gmycFactors, ".log", sep=""))
## ## stopifnot(all(file.exists(logFiles)))
## ## resTmrca <- lapply(logFiles, function(logF) {
## ##     tmpDt <- read.table(logF, header=TRUE, colClasses="numeric")
## ##     stopifnot(nrow(tmpDt) == 10001)
## ##     colid <- grep("tmrca", names(tmpDt))
## ##     stopifnot(length(colid) == 1)
## ##     rg <- 1001:nrow(tmpDt)
## ##     list(mean=mean(tmpDt[rg, colid]),
## ##          median=median(tmpDt[rg, colid]))
## ## })
## ## names(resTmrca) <- logFiles
## ## save(resTmrca, file="data/resTmrca.RData")
## load(file="data/resTmrca.RData")
## meanTmrca <- sapply(resTmrca, function(x) x$mean)
## medianTmrca <- sapply(resTmrca, function(x) x$median)

## ### starBEAST PPS
## ### ---- starbeast-pps ----
## library(starbeastPPS)

## impxml <- read.starbeast(beast.xml="/home/francois/Photos/impatiens_analyses/000.allESU1/starbeastPPS_allESU1_longSampling_run1/20131230.allESU1_longSampling.xml", "/home/francois/Photos/impatiens_analyses/000.allESU1/starbeastPPS_allESU1_longSampling_run1/combined")

## impcoal <- analyze.coalescent(impxml, "~/Software/msdir/")

## impseq <- analyze.sequences(impxml, "/home/francois/Software/Seq-Gen.v1.3.3")

## ### bGMYC analysis
## ### ---- bgmyc ----
## library(bGMYC)

## resSingle <- bgmyc.singlephy(trBeast, mcmc=50000, burnin=1, thinning=10, t1=2, t2=100, start=c(1,1,25))

## ### ---- ClockStaR ----
## ## system("cd ~/R-dev/; git clone git@github.com:sebastianduchene/ClockstaR.git")
## ## install("~/R-dev/ClockstaR")
## library(ClockstaR2)
## optim.trees.interactive()

## #####################

## ### Just for ESU1
## ## impESU1Alg <- mergeSeq(impESU1, output="/tmp", seqFolder="~/Documents/seqRepository",
## ##                        markers=c("16Sc", "16S", "COI", "ATP6", "H3a", "LSU"))
## ## concatenateAlignments(patter="afa$", path="/tmp", output="/tmp/20130202.impESU1.phy",
## ##                       partition="/tmp/20130202.impESU1.part", partition.format="nexus",
## ##                       create.conc=T, colsep="", colw=10000)


## ### to get all extractions that have at least another marker other than COI
## ## impAlg <- mergeSeq(impExt, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
## ##                    markers=c("16S", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a", "COI"),
## ##                    justCheck=T)
## ## impAlg <- impAlg[(impAlg$"16S" | impAlg$ATP6 | impAlg$c0036 | impAlg$ITS | impAlg$LSU | impAlg$H3a), ]
## ## impExt2 <- rownames(impAlg)
## ## impAlg <-  mergeSeq(impExt2, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
## ##                    markers=c("16S", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a", "COI"))
## ## concatenateAlignments(patter="afa$", path="/tmp/seq", output="/tmp/20130411.impatiens.phy",
## ##                       partition="/tmp/20130411.impatiens.part", partition.format="nexus",
## ##                       create.conc=T, colsep="", colw=10000)

## ## ### compare individual topologies (NJ based)
## ## trc0036 <- read.tree(file="20130507-c0036.tre")
## ## trc0775 <- read.tree(file="20130507-c0775.tre")
## ## trLSU <- read.tree(file="20130507-LSU.tre")
## ## trITS <- read.tree(file="20130507-ITS.tre")
## ## trH3a <- read.tree(file="../20130507-H3a.tre")

## ## trc0036 <- extract_to_label(trc0036, impDB)
## ## trc0775 <- extract_to_label(trc0775, impDB)
## ## trLSU <- extract_to_label(trLSU, impDB)
## ## trITS <- extract_to_label(trITS, impDB)
## ## trH3a <- extract_to_label(trH3a, impDB)

## ## pdf(file="/tmp/trees.pdf")
## ## plot(trc0036, no.margin=T, cex=.5)
## ## plot(trc0775, no.margin=T, cex=.5)
## ## plot(trLSU, no.margin=T, cex=.5)
## ## plot(trITS, no.margin=T, cex=.5)
## ## plot(H3a, no.margin=T, cex=.5)
## ## dev.off()


## ### run RAxML to get reduced file
## ## phy2nex(file="~/Documents/Impatiens/20130507.impatiens_partfinder/20130507.impatiens.phy",
## ##         partition.file="~/Documents/Impatiens/20130507.impatiens_partfinder/20130507.impatiens.part")

## ## rownames(impAlg) <- extract_to_labelStr(rownames(impAlg), impDB)

## ## tmpTr <- read.tree(file="/tmp/seq/RAxML_bestTree.T1")
## ## tmpTr <- root(tmpTr, grep("S0213", tmpTr$tip.label))
## ## tmpTr$tip.label <- extract_to_labelStr(tmpTr$tip.label, impDB)
## ## tmpTr <- ladderize(tmpTr)
## ## posMarkers <- seq(1, by=.02, length.out=ncol(impAlg))
## ## tipOrder <- tmpTr$tip.label[tmpTr$edge[which(tmpTr$edge %in% 1:Ntip(tmpTr))]]
## ## pdf(file="/tmp/20121022.tree.pdf", height=20)
## ## plot(tmpTr, no.margin=TRUE, cex=.4)
## ## points(rep(posMarkers[1], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder, ]$"16Sc", 15, 22), cex=.6)
## ## points(rep(posMarkers[2], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"16S", 15, 22), cex=.6)
## ## points(rep(posMarkers[3], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"COI", 15, 22), cex=.6)
## ## points(rep(posMarkers[4], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"ATP6", 15, 22), cex=.6)
## ## points(rep(posMarkers[5], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"H3a", 15, 22), cex=.6)
## ## points(rep(posMarkers[6], Ntip(tmpTr)), 1:Ntip(tmpTr), pch=ifelse(impAlg[tipOrder,]$"LSU", 15, 22), cex=.6)
## ## text(posMarkers, rep(Ntip(tmpTr)+1, length(posMarkers)), c("16Sc", "16S", "COI", "ATP6", "H3a", "LSU"), srt=90, cex=.4)
## ## dev.off()


## ### summary table
## ## impAlgNm <- extract_to_string_label(rownames(impAlg), impDB)
## ## impAlg$ESU <- sapply(impAlgNm, function(x) unlist(strsplit(x, "_"))[3])

## ## library(reshape)
## ## impSum <- melt(impAlg, id.vars="ESU", variable_name="marker")

## ## xtabs(~ ESU + marker, data=impSum, subset=value)

## ### haplotype networks
##  library(pegas)
## esu1HaploData <- subset(impDB, consensusESU %in% c("ESU1"))$Extract
## esu1HaploData <- esu1HaploData[nzchar(esu1HaploData)]
## esu1HaploData <- gsub(",.+$", "", esu1HaploData)

## esu2HaploData <- subset(impDB, consensusESU %in% c("ESU2"))$Extract
## esu2HaploData <- esu2HaploData[nzchar(esu2HaploData)]
## esu2HaploData <- gsub(",.+$", "", esu2HaploData)

## impCOI <- read.dna(file="data/20130923.impatiens_COI.phy", format="seq")

## esu1Seq <- impCOI[match(esu1HaploData, dimnames(impCOI)[[1]]), ]
## esu1ToRm <- which(sapply(esu1Seq, function(x) any(is.na(base.freq(x)))))
## esu1Seq <- esu1Seq[-esu1ToRm, ]

## esu2Seq <- impCOI[match(esu2HaploData, dimnames(impCOI)[[1]]), ]
## esu2ToRm <- which(sapply(esu2Seq, function(x) any(is.na(base.freq(x)))))
## esu2Seq <- esu2Seq[-esu2ToRm, ]


## esu1Haplotypes <- haplotype(esu1Seq)
## esu1HaploNet <- haploNet(esu1Haplotypes)
## plot(esu1HaploNet, size=attr(esu1HaploNet, "freq"), bg=ziss, scale.ratio=20, cex=.6)


## esu2Haplotypes <- haplotype(esu2Seq)
## esu2HaploNet <- haploNet(esu2Haplotypes)
## plot(esu2HaploNet, size=attr(esu2HaploNet, "freq"))

## ### test GMYC
## ## library(splits)
## ## trBeast <- read.nexus(file="/home/francois/Documents/Impatiens/20130201.beast_COI_PB_strictClock/20130201-COI-PB-strict.tre")
## ## singleGmycBeast <- gmyc(trBeast)
## ## multiGmycBeast <- gmyc(trBeast, method="multiple")

## ## trNJ <- read.tree(file="~/Documents/Impatiens/20121105.impatiens_beast/20130131.COI-NJ_ultra.nex")
## ## trNJ <- drop.tip(trNJ, "S0213")
## ## trNJ <- multi2di(trNJ)
## ## singleGmycNJ <- gmyc(trNJ)
## ## multiGmycNJ <- gmyc(trNJ, method="multiple")

## ## trBeastCOI <- read.nexus(file="~/Documents/Impatiens/20130131.impatiens_beast_COI/20130131-111957-COI_tree.nex")
## ## trBeastCOI <- drop.tip(trBeastCOI, "S0213")
## ## singleGmycBeastCOI <- gmyc(trBeastCOI)
## ## multiGmycBeastCOI <- gmyc(trBeastCOI, method="multiple")

## #####################
