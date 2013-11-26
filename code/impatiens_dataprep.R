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

makeTraitFile <- function(alg, output, db=impDB, allowedESUs,
                          format="fasta", ...) {
### generate the traits file to be used by *BEAST
### alg -- path to alignment in PHYLIP/FASTA format
### output -- path and filename of where to write the output
### allowedESUs -- safe guard to validate that the ESU in database match
####   what is intended.
### ... -- additional arguments to be passed to ape:::read.dna
  if (!file.exists(alg)) stop(alg, "doesn't exist")
  algF <- read.dna(file=alg, format=format, ...)
  tmpRes <- cbind(traits=db$Extract, species=db$consensusESU)
  tmpRes <- tmpRes[nzchar(tmpRes[,1]), ]
  tmpRes[,1] <- gsub(",.+$", "", tmpRes[,1]) # replace extractions by first one if multiple to follow mergeSeq rule
  tmpRes <- subset(tmpRes, tmpRes[,1] %in% dimnames(algF)[[1]])
  notAllowed <- tmpRes[,2][! tmpRes[,2] %in% allowedESUs]
  if (length(notAllowed) > 0) warning(paste(unique(notAllowed), collapse=","), ", not found in DB")
  if (any(is.na(tmpRes[,2]))) warning("some NAs found in DB")
  if (any(!nzchar(tmpRes[,2]))) warning("some ESUs not specified")
  write.table(tmpRes, file=output, row.names=FALSE, quote=FALSE, sep="\t")
  TRUE
}


impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv", stringsAsFactors=FALSE)
impExt <- impDB$Extract[nzchar(impDB$Extract)]
stopifnot(ncol(impDB) > 1)

### Main analysis -- these commands should produce the final output for the published
###  analyses. If it's not the case, it should go somewhere else.

## Create alignments for individual markers for main phylogeny
impAlg <- mergeSeq(impExt, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
                   markers=c("16S", "16Sc", "COI", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a"),
                   gblocks=list("ITS" = "-t=d -b4=5 -b5=a -p=n"), justCheck=F)
concatenateAlignments(pattern="afa$", path="/tmp/seq", output="/tmp/seq/20130923.impatiens.phy",
                      partition="/tmp/seq/20130923.impatiens.part", partition.format="nexus",
                      create.conc=TRUE, colsep="", colw=10000)
alg2nex(file="/tmp/seq/20130923.impatiens.phy", partition.file="/tmp/seq/20130923.impatiens.part")

## unfortunately at this stage, still need manual editing of the partition... compare -orig.nex and .nex

## Create alignements for individual markers for *BEAST analysis of ESU1 "complex"
##   For this, first look at available sequences for the complex and select specimens
##   that have at least 1 mt marker and 1 nuc marker (except for LSU that is too conserved).

impExtESU1tmp <- subset(impDB, consensusESU %in% c("ESU1", "ESU3", "gracilis", "Hawaii", "Wpac", "RedSea"))
impExtESU1tmp <- impExtESU1tmp$Extract[nzchar(impExtESU1tmp$Extract)]
impAlgESU1tmp <- mergeSeq(impExtESU1tmp, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
                          markers=c("16S", "16Sc", "COI", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a"),
                          justCheck=TRUE)
impExtESU1 <- rownames(subset(impAlgESU1tmp, (impAlgESU1tmp$"16S" == TRUE | impAlgESU1tmp$"16Sc" == TRUE |
                                              impAlgESU1tmp$"COI" == TRUE | impAlgESU1tmp$"ATP6" == TRUE) &
                              (impAlgESU1tmp$"c0036" == TRUE | impAlgESU1tmp$"c0775" == TRUE |
                               impAlgESU1tmp$"ITS" == TRUE | impAlgESU1tmp$"H3a" == TRUE)))
impAlgESU1 <- mergeSeq(impExtESU1, output="/tmp/seq", seqFolder="~/Documents/seqRepository",
                       markers=c("16S", "16Sc", "COI", "ATP6", "c0036", "c0775", "ITS", "LSU", "H3a"),
                       gblocks=list("ITS" = "-t=d -b4=5 -b5=a -p=n"), justCheck=F)
concatenateAlignments(pattern="afa$", path="~/Documents/Impatiens/20130929.impatiens_sampling",
                      output="~/Documents/Impatiens/20130929.impatiens_sampling/20130930.impatiens_starbeast.phy",
                      partition="~/Documents/Impatiens/20130929.impatiens_sampling/partition",
                      partition.format="nexus", colw=10000, colsep="")
alg2nex(file="~/Documents/Impatiens/20130929.impatiens_sampling/20130930.impatiens_starbeast.phy",
        partition.file="~/Documents/Impatiens/20130929.impatiens_sampling/partition")

makeTraitFile(alg="/tmp/seq/20130930-133004-COI.afa", output="/tmp/seq/traitsAllESUs.txt",
              allowedESUs=c("ESU1", "ESU3", "RedSea", "Hawaii", "Wpac", "gracilis"))

## Redo analyses without mtDNA
concatenateAlignments(pattern="(c0036|c0775|H3a|ITS|LSU).+afa$",
                      path="~/Documents/Impatiens/20130929.impatiens_sampling/origAlignments",
                      output="~/Documents/Impatiens/20130929.impatiens_sampling/20130929.impatiens_noNuc_allESU1.phy",
                      partition="/tmp/partitionnuc", partition.format="nexus", create.conc=TRUE, colsep="", colw=10000)

alg2nex(file="~/Documents/Impatiens/20130929.impatiens_sampling/20130929.impatiens_noNuc_allESU1.phy",
        partition.file="/tmp/partitionnuc")

## fasToPhase("/tmp/seq/20130929-172701-c0036.afa")
## fasToPhase("/tmp/seq/20130929-172701-c0775.afa")
## fasToPhase("/tmp/seq/20130929-172701-ITS.afa")

## system("PHASE -d1 /tmp/seq/20130929-172701-c0036.inp /tmp/seq/20130929-172701-c0036.out")
## system("PHASE -d1 /tmp/seq/20130929-172701-c0775.inp /tmp/seq/20130929-172701-c0775.out")
## system("PHASE -d1 /tmp/seq/20130929-172701-ITS.inp /tmp/seq/20130929-172701-ITS.out")

## phaseToFas("/tmp/seq/20130929-172701-c0036.out", "/tmp/seq/20130929-172701-c0036.afa",
##            "/tmp/seq/20130929-172701-c0036-phased.afa", removeEmptySeqs=TRUE)

## phaseToFas("/tmp/seq/20130929-172701-c0775.out", "/tmp/seq/20130929-172701-c0775.afa",
##            "/tmp/seq/20130929-172701-c0775-phased.afa", removeEmptySeqs=TRUE)

## phaseToFas("/tmp/seq/20130929-172701-ITS.out", "/tmp/seq/20130929-172701-ITS.afa",
##            "/tmp/seq/20130929-172701-ITS-phased.afa", removeEmptySeqs=TRUE) ## problematic output; only 1 ambiguity
