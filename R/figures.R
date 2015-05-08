tikz_ <- function(obj, file, width, height, standAlone = TRUE, ...) {
    tikzDevice::tikz(file = file, width = width, height = height,
                     packages = c("\n\\nonstopmode\n", getOption("tikzXelatexPackages")),
                     standAlone = standAlone, ...)

    on.exit(dev.off())
    if (inherits(obj, "ggplot"))
        print(obj)
    else {
        eval.parent(substitute(obj))
    }
}


plot2pdf <- function(plt, pdf.output, width, height, latex = "xetex", ...) {
    stopifnot(grepl("\\.pdf$", pdf.output))
    unlink(tikz.plot <- gsub("\\.pdf$", ".tikz", pdf.output))
    unlink(log.plot <- gsub("\\.pdf$", ".log", pdf.output))

    tikz_(plt, file = tikz.plot, width = width, height = height, ...)

    ## from knitr
    unlink(pdf.output)
    latex <- switch(latex,
                    pdftex = getOption('tikzLatex'),
                    xetex  = getOption('tikzXelatex'),
                    luatex = getOption('tikzLualatex'),
                    stop('a LaTeX engine must be specified for tikzDevice', call. = FALSE)
                    )
    owd <- setwd(dirname(pdf.output))
    system2(latex, shQuote(basename(tikz.plot)), stdout = NULL)
    setwd(owd)

    if (!file.exists(pdf.output)) {
        if (file.exists(log.plot)) {
            message(paste(readLines(log.plot), collapse = "\n"))
        }
        stop("failed to compile ", tikz.plot, ' to PDF',  call. = FALSE)
    }
    return(NULL)
}

gmyc_coi_plot <- function(gmyc_summary, ...) {
    levels(gmyc_summary$sequences)[levels(gmyc_summary$sequences) == "allSeq"] <- "All Sequences"
    levels(gmyc_summary$sequences)[levels(gmyc_summary$sequences) == "noDup"]  <- "Unique haplotypes"
    levels(gmyc_summary$clock)[levels(gmyc_summary$clock) == "strict"] <- "Strict clock"
    levels(gmyc_summary$clock)[levels(gmyc_summary$clock) == "relaxed"] <- "Relaxed clock"
    p <- ggplot(data=gmyc_summary, aes(x=demographic, y=mean,
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
      scale_color_manual(values = wesanderson::wes_palette("Zissou", 5)[c(1, 5)],
                         labels=c("multi-threshold GMYC", "single threshold GMYC"))
    plot2pdf(p, ...)
}


starbeast_summary_plot <- function(star_beast_summary_all, ...) {
    sbSummPlot <- reshape2::melt(star_beast_summary_all[, -match(c("PS_logLik", "SS_logLik", "stdSS", "stdPS"),
                                                                 names(star_beast_summary_all))],
                                 measure.vars=c("BFSS", "BFPS"))

    sbCol <- wesanderson::wes_palette("Zissou", 5)[c(1, 5)]
    names(sbCol) <- c("BFSS", "BFPS")

    p <- ggplot(sbSummPlot, aes(x=plotGroupings, y=value, group=variable, colour=variable)) +
      scale_x_discrete(limits=c("M0", "oversplit", "M5", "M4", "M6", "M3", "M2", "M1")) +
      scale_colour_manual(values=sbCol, breaks=c("BFSS", "BFPS"),
                          labels=c("Stepping Stone Sampling", "Path Sampling")) +
      geom_point(position=position_dodge(width=.3)) +
      stat_summary(fun.y = mean, geom="point", size=3,
                   position=position_dodge(width=.3)) +
      geom_hline(yintercept=-10, width=.1, col=sbCol[2], linetype=2) +
      theme(legend.position="top", legend.title=element_blank(),
            panel.background = element_rect(fill = "gray95"),
            panel.grid.major = element_line(colour = "white", size=0.1),
            panel.grid.minor = element_line(NA)) +
      labs(y="Bayes Factors ", x="Models") +
      theme(legend.justification=c(0,0), legend.position=c(0,0)) +
      facet_grid( ~ dataIncluded)

    plot2pdf(p, ...)
}

loci_coverage_plot <- function(locus_graph_data, ...) {
    lbl <- c("16S", "COI", "ATP6", "c0036", "c0775", "H3a", "ITS", "LSU")
    sub_loc <- locus_graph_data[match(lbl, locus_graph_data[["Locus"]]), ]
    sub_loc$Center <- (sub_loc$end + sub_loc$begin)/2
    p <- ggplot(data=locus_graph_data) +
      geom_segment(aes(x=begin, xend=end, y=Order, yend=Order, colour=consensusESU),
                   lineend="round", size=I(1.5)) +
      scale_x_continuous(breaks=sub_loc$Center, labels=lbl) +
      scale_y_discrete(labels=element_blank()) + ylab("Individuals") + xlab("Loci") +
      scale_colour_manual(values=load_impPal()) +
      theme(legend.position=c(0, 1),
            panel.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.ticks=element_blank())
    plot2pdf(p, ...)
}


impatiens_tree_plot <- function(impTree, impDB, ...) {
    imp_tree_plot_ <- function(impTree, impDB, impPal) {
        impTree <- ladderize(impTree)
        posTips <- max(branching.times(impTree))
        impTree <- extract_to_label(impTree, impDB, c("consensusESU", "Country", "UFID", "Extract"))
        esuList <- c("ESU1", "ESU2", "ESU3", "gracilis", "tiger", "tigerRedSea", "Medit", "WA",
                     "Gala", "EP", "Hawaii", "Wpac", "RedSea")

        impNodLbl <- impTree$posterior
        impNodLblCol <- rep(NULL, length(impNodLbl))
        impNodLblCol[impNodLbl >= .99] <- "black"

        impAllNds <- 1:length(impNodLbl)
        impKeepNds <- impAllNds[!is.na(impNodLblCol)]

        par(mai=c(1,0,0,0))
        plot.phylo(impTree, root.edge=TRUE, show.tip.label=FALSE, x.lim=c(0, 12.5))
        nodelabels(text=rep("", length(impKeepNds)), node=impKeepNds+Ntip(impTree),
                   frame="circ", col=impNodLblCol[!is.na(impNodLblCol)],
               bg=impNodLblCol[!is.na(impNodLblCol)],
                   ##fg=impNodLblCol[!is.na(impNodLblCol)],
                   cex=.3)
        barMonophyletic(groupLabel=esuList, groupMatch=paste("^", esuList, "_", sep=""), impTree, cex.text=.4,
                        cex.plot=.3, extra.space=.1, text.offset=1.02,
                        seg.col=load_impPal()[esuList])
        abline(v=max(branching.times(impTree)) - c(1:4, 6:10), lty=1, col="gray70")
        abline(v=max(branching.times(impTree) - 5), lty=2, col="gray50")
        axis(side=1, at=max(branching.times(impTree)) - 0:11, labels=c(0, paste("-", 1:10, sep=""), NA))
        legend(x=0, y=10, pch=16, col="black", legend=c("PP $\\geq$ 0.99"), bg="white")
    }
    plot2pdf(imp_tree_plot_(impTree = impTree, impDB), ... )
}

## 1. WA + EP + Gala
draw_map_WA <- function(impatiens_map_data, ...) {
    p <- draw_impatiens_map(impatiens_map_data, bg = map_data("world2"),
                            palette = load_impPal(), ESUs =  c("WA", "EP", "Gala"),
                            x.lim = c(250, 310), y.lim = c(-25, 25))
    plot2pdf(p, ...)
}


## 2. tiger + ESU2 + redSeaTiger
draw_map_group2 <- function(impatiens_map_data, ...) {
    p <- draw_impatiens_map(impatiens_map_data, bg = map_data("world2"),
                            palette = load_impPal(),
                            ESUs =  c("tiger", "ESU2", "tigerRedSea"),
                            x.lim =c(25, 220), y.lim = c(-25, 25))
    plot2pdf(p, ...)
}


## 3. ESU3 + RedSea + Gracilis + Hawaii + WPac + ESU1
draw_map_group1 <- function(impatiens_map_data, ...) {
    p <- draw_impatiens_map(impatiens_map_data, bg = map_data("world2"),
                            palette = load_impPal(),
                            ESUs = c("ESU1", "ESU3", "RedSea", "gracilis", "Hawaii", "Wpac"),
                            x.lim = c(30, 220), y.lim = c(-27, 29))
    plot2pdf(p, ...)
}


raxml_tree_plot <- function(raxml_tree_file, impDB, ...) {

    raxml_tree_ <- function(file, impDB) {
        raxmlImpTr <- ape::read.tree(file=file)
        raxmlImpTr <- ape::drop.tip(root(raxmlImpTr, "S0213"), "S0213")
        raxmlImpTr <- ape::ladderize(raxmlImpTr)

        raxmlImpTr <- extract_to_label(raxmlImpTr, impDB, c("consensusESU", "Country", "UFID", "Extract"))

        ndCol <- character(length(raxmlImpTr$node.label))
        ndCol[as.numeric(raxmlImpTr$node.label) >= 90] <- "black"
        ndCol[as.numeric(raxmlImpTr$node.label) >= 80 &
                as.numeric(raxmlImpTr$node.label) < 90] <- "red"
        ndAll <- 1:length(ndCol)
        ndKeep <- ndAll[nzchar(ndCol)]

        par(mai=c(0,0,0,0))
        ape::plot.phylo(raxmlImpTr, cex=.6, show.tip.label=FALSE, x.lim=c(0,.245))
        chopper::barMonophyletic(groupLabel=load_esuList(),
                                 groupMatch=paste("^", load_esuList(), "_", sep=""),
                                 raxmlImpTr, cex.text=.4, cex.plot=.3, extra.space=.155,
                                 text.offset=1.02, seg.col=load_impPal()[load_esuList()])
        ape::nodelabels(text=rep("", length(ndKeep)), node=ndKeep+Ntip(raxmlImpTr),
                        bg=ndCol[nzchar(ndCol)], frame="circ", cex=.275)
        ape::add.scale.bar()
    }
    plot2pdf(raxml_tree_(file = raxml_tree_file, impDB), ...)
}

per_locus_trees <- function(raxFiles, impDB, ...) {
    uniqExt <- regmatches(impDB$Extract, regexpr("^[^,]+", impDB$Extract))
    loci <- c("c0036", "c0775", "H3a", "mtDNA", "rDNA")

    per_locus_trees_ <- function(raxFiles, impDB) {
        par(mai=c(0,0,.5,0))
        layout(matrix(c(4,1,2,4,3,5), 2, 3, byrow=TRUE))
        for (i in 1:length(raxFiles)) {
            raxmlTr <- read.tree(file=raxFiles[i])
            raxmlCol <- impDB[match(raxmlTr$tip.label, uniqExt), "consensusESU"]
            raxmlCol <- load_impPal()[raxmlCol]
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
                legend(x=0, y=70, legend=textLgd, col=load_impPal()[textLgd], lty=1, lwd=3)
            }
        }
    }
    plot2pdf(per_locus_trees_(raxFiles, impDB), ...)

}

get_gmyc_results_names <- function(gmyc_results) {
    nmTrees <- gsub(".+[0-9]\\.(.+)\\..+\\..+$", "\\1", names(gmyc_results))
    nmTrees <- gsub("^allSeq_", "All sequences, ", nmTrees)
    nmTrees <- gsub("^noDup_", "Haplotypes, ", nmTrees)
    nmTrees <- gsub("relaxed_", "Relaxed Clock, ", nmTrees)
    nmTrees <- gsub("strict_", "Strict Clock, ", nmTrees)
    nmTrees <- gsub("yule$", "Yule", nmTrees)
    nmTrees <- gsub("coalexp", "Coalescent (Exponential growth)", nmTrees)
    nmTrees <- gsub("coalcst", "Coalescent (Constant)", nmTrees)
    nmTrees
}

gmyc_tree_all_sequences <- function(gmyc_results, ...) {

    gmyc_tree_all_ <- function(gmyc_results) {
        nmTrees <- get_gmyc_results_names(gmyc_results)
        allSeqGmEd <- lapply(gmyc_results[1:6], function(x) gmycEdges(x$simpleGmyc))

        par(mfrow=c(3, 2))
        for (i in 1:6) {
            resKeep <- gmyc_results[[i]]$simpleGmyc
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
                colToUse <- load_impPal()[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3", "gracilis",
                                            "tiger", "Gala", "Medit", "ESU3_Deep", "ESU1_Lizard", "RedSea",
                                            "Wpac")]
            } else if (nSpp == 16) {
                colToUse <- load_impPal()[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3",
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
    }
    plot2pdf(gmyc_tree_all_(gmyc_results), ...)
}

gmyc_tree_unique_sequences <- function(gmyc_results, ...) {

    gmyc_tree_uniq_ <- function(gmyc_results) {
        nmTrees <- get_gmyc_results_names(gmyc_results)
        uniqSeqGmEd <- lapply(gmyc_results[7:12], function(x) gmycEdges(x$simpleGmyc))

        par(mfrow=c(3, 2))
        for (i in 7:12) {
            k <- i - 6
            resKeep <- gmyc_results[[i]]$simpleGmyc
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
                colToUse <- load_impPal()[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3", "gracilis",
                                            "tiger", "Gala", "Medit", "ESU3_Deep", "ESU1_Lizard", "RedSea",
                                            "Wpac")]
            }  else if (nSpp == 16) {
                colToUse <- load_impPal()[c("EP", "Hawaii", "WA", "ESU1", "ESU2", "ESU3",
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
    }
    plot2pdf(gmyc_tree_uniq_(gmyc_results), ...)
}
