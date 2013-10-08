
##' Reformat the output from BEAST to make it readable by Tracer
##'
##' My runs got messed up as I asked to log way too many states. I
##' ended up with corrupted, > 2Gb log files. To fix the problem, I'm
##' parsing out the output of BEAST, saved in STDOUT, and recreate a
##' log file and a file that contains the marginal likelihoods from
##' path sampling and stepping stone sampling.
##' @title Reformat BEAST stdout
##' @param file Full path for the file that contains de stdout
##' @param outputLog Full path for the newly created log file
##' @param outputMarginal Full path for the file containing the
##' marginal likelihoods
##' @return nothing, just creates files.
##' @author François Michonneau
reformatBeastStdOut <- function(file, outputLog, outputMarginal) {

    stdout <- scan(what="character", sep="\n", file=file)

    if (!is.null(outputLog)) {
        hasNum <- sapply(stdout, function(x) length(grep("^[0-9]", x)) > 0)
        header <- which(hasNum)[1] - 1
        cat(stdout[header], stdout[hasNum], file=outputLog, sep="\n")
    }

    if (!is.null(outputMarginal)) {
        logLik <- sapply(stdout, function(x) length(grep("marginal", x)) > 0)
        cat(stdout[logLik], file=outputMarginal, sep="\n")
    }
}

reformatBeastStdOut(file="~/Photos/impatiens_analyses/20131002.allImpatiens/stdout.txt",
                    outputLog="~/Photos/impatiens_analyses/20131002.allImpatiens/newLog.txt",
                    outputMarginal="~/Photos/impatiens_analyses/20131002.allImpatiens/marginalLik.txt")

reformatBeastStdOut(file="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaii/stdout.txt",
                    outputLog="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaii/newLog.txt",
                    outputMarginal="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaii/marginalLik.txt")

reformatBeastStdOut(file="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaiiRedSea/stdout.txt",
                    outputLog="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaiiRedSea/newLog.txt",
                    outputMarginal="~/Photos/impatiens_analyses/20131002.allImpatiens_noWpacHawaiiRedSea/marginalLik.txt")

reformatBeastStdOut(file="~/Photos/impatiens_analyses/20131002.allImpatiens_noHawaii/impbeast.out",
                    outputLog=NULL,
                    outputMarginal="~/Photos/impatiens_analyses/20131002.allImpatiens_noHawaii/marginalLik.txt")
