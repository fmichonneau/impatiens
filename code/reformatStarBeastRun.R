
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


reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.allESU1/20131007.impatiens_allESU1/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.allESU1/20131007.impatiens_allESU1/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noHawaii/20131007.impatiens_noHawaii/impbeast.out",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noHawaii/20131007.impatiens_noHawaii/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noRedSea/20131002.impatiens_noRedSea_run1/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noRedSea/20131002.impatiens_noRedSea_run1/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noRedSea/20131007.impatiens_noRedSea_run2/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noRedSea/20131007.impatiens_noRedSea_run2/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noWpac/20131007.impatiens_noWpac/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noWpac/20131007.impatiens_noWpac/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noWpacHawaii/20131002.impatiens_noWpacHawaii/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noWpacHawaii/20131002.impatiens_noWpacHawaii/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noWpacHawaii/20131007.impatiens_noWpacHawaii/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noWpacHawaii/20131007.impatiens_noWpacHawaii/marginalLik.txt")


reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noWpacHawaiiRedSea/20131007.impatiens_noWpacHawaiiRedSea/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noWpacHawaiiRedSea/20131007.impatiens_noWpacHawaiiRedSea/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.noWpacRedSea/20131007.impatiens_noWpacRedSea/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.noWpacRedSea/20131007.impatiens_noWpacRedSea/marginalLik.txt")

reformatBeastStdOut(file="/home/francois/Documents/Impatiens/000.random/20131007.impatiens_random/stdout.txt",
                    outputLog=NULL,
                    outputMarginal="/home/francois/Documents/Impatiens/000.random/20131007.impatiens_random/marginalLik.txt")


##################


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
