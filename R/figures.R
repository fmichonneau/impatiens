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

