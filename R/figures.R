plot2pdf <- function(plt, pdf.output, width, height, latex = "xetex", ...) {
    stopifnot(grepl("\\.pdf$", pdf.output))
    unlink(tikz.plot <- gsub("\\.pdf$", ".tikz", pdf.output))
    unlink(log.plot <- gsub("\\.pdf$", ".log", pdf.output))

    tikz(file = tikz.plot, width = width, height = height,
         packages = c("\n\\nonstopmode\n", getOption("tikzXelatexPackages")), ...)
    print(plt)
    dev.off()

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
