get_locus_files <- function(algfile, partfile) {
    locFiles <- chopper::cutAlignment(algfile = algfile, partfile = partfile,
                                      formatin="sequential", format="sequential")
    part_info <- chopper::raxmlPartitionInfo(partfile)
    attr(locFiles, "partition_info") <- part_info
    locFiles
}

extract_locus_characteristics <- function(locus_file, tmpDir) {
    tmpfnm <- file.path(tmpDir, basename(locus_file))
    chopper::convertGaps(locus_file, output=tmpfnm, formatin="sequential", colw=10000)
    chopper::removeEmptySeqs(tmpfnm, formatin="sequential", formatout="sequential",
                             gap="-", overwrite=TRUE, colw=10000)
    algChar <- ape::read.dna(file=tmpfnm, format="sequential", as.character=TRUE)
    alg <- ape::read.dna(file=tmpfnm, format="sequential")
    lSeqNonAlign <- apply(algChar, 1, function(x) {
        seq <- paste(as.character(x), collapse="", sep="")
        length(gregexpr("[actgACTG]", seq)[[1]])
    })
    nSeq <- nrow(alg)                    # number of individuals sequenced
    nUniq <- nrow(haplotype(alg))        # number of haplotypes/unique sequences
    lSeq <- ncol(alg)                    # length of the aligned sequence
    rglSeq <- range(lSeqNonAlign)        # range of raw sequence length
    nSeg <- length(ape::seg.sites(alg))  # number of segregating/variable sites
    nPis <- phyloch::pis(alg)            # number of parsimony informative sites
    c(nSeq, nUniq, lSeq, rglSeq[2], nSeg, nPis)
}

extract_loci_characteristics <- function(locFiles, tmpDir = "tmp") {

    if (! file.exists(tmpDir)) {
            dir.create(tmpDir)
    }

    lociChar <- sapply(locFiles, extract_locus_characteristics, tmpDir)
    lociChar <- as.matrix(lociChar)
    locNm <- gsub(".+_(.+)\\.phy", "\\1", colnames(lociChar))
    colnames(lociChar) <- locNm
    rownames(lociChar) <- c("Nind",
                            "K",
                            "bp_alg",
                            "bp_unalg",
                            "S", "S_i")
    lociChar
}

generate_loci_table <- function(lociChar) {

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
                    "\\\\multicolumn{5}{|c}{nucDNA} \\\\\\\\ \\\n",
                    "\\\\cline{2-4} \\\\cline{5-9}", sep="")

    cat(gsub("(}\\\n)(\\s+\\\\\\hline)", multiColStr, locTable))
}

get_locus_coverage <- function(locus_files) {
    locCov <- sapply(locus_files, function(fnm) {
                         tmpfnm <- file.path("tmp", basename(fnm))
                         chopper::convertGaps(fnm, output=tmpfnm, formatin="sequential", colw=10000)
                         chopper::removeEmptySeqs(tmpfnm, formatin="sequential", formatout="sequential",
                                                  gap="-", overwrite=TRUE, colw=10000)
                         alg <- ape::read.dna(file=tmpfnm, format="sequential")
                         dimnames(alg)[[1]]
                     })
    names(locCov) <- gsub(".+_(.+)\\.phy", "\\1", names(locCov))
    locCov
}

calculate_locus_positions <- function(loc_info, pos) {
    res <- vector("list", length(pos))
    for (i in seq_along(pos)) {
        if (i == 1) {
            beg <- 1
        } else {
            beg <- max(unlist(res)) + 1
        }
        end <- (beg + loc_info[pos[i]]) - 1
        res[[i]] <- c(begin = beg, end = end)
    }
    names(res) <- names(loc_info)[pos]
    res
}


get_locus_graph_data <- function(locus_coverage, impDB, locus_files) {
    locWhich <- mapply(function(x, y) cbind(x, rep(y, length(x))),
                       locus_coverage, names(locus_coverage))
    locWhich <- do.call("rbind", locWhich)
    locWhich <- data.frame(locWhich, stringsAsFactors=FALSE)
    names(locWhich) <- c("Extract", "Locus")

    part_info <- attr(locus_files, "partition_info")
    locPos <- sapply(part_info, function(x) x[[1]][["end"]] - x[[1]][["begin"]])
    locPos <- calculate_locus_positions(locPos, pos = c(1, 5, 2, 3, 4, 6, 7, 8))
    locPos <- do.call("rbind", locPos)
    locPos <- as.data.frame(locPos)
    names(locPos) <- c("begin", "end")
    locPos$Locus <- rownames(locPos)

    locWhich <- merge(locWhich, locPos)

    tmpDB <- impDB[, c("Extract", "consensusESU")]
    tmpDB <- tmpDB[nzchar(tmpDB$Extract), ]
    tmpDB$Extract <- sapply(tmpDB$Extract, function(x) unlist(strsplit(x, ","))[1])
    locWhich <- merge(locWhich, tmpDB)

    locWhich <- subset(locWhich, Extract != "S0213")

    locRank <- data.frame(table(locWhich$Extract))
    names(locRank) <- c("Extract", "nLoci")
    locRank$Rank <- rank(locRank$nLoci, ties.method="random")
    locWhich <- merge(locWhich, locRank)

    locWhich <- locWhich[order(locWhich$nLoci, locWhich$consensusESU), ]
    ordr <- data.frame(Extract=unique(locWhich$Extract), Order=1:length(unique(locWhich$Extract)))
    locWhich <- merge(locWhich, ordr)
    locWhich
}
