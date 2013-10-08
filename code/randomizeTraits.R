##' Randomize a *BEAST trait file
##'
##' This function takes a *BEAST trait file and randomize the species
##' assignment. It outputs the relevant XML parts that need to be
##' replaced in the original file.
##' @title Randomize *BEAST trait file
##' @param traitFile a character string indicating the path of the
##' original *BEAST trait file.
##' @param output a character string indicating the path of the output file.
##' @return Nothing. Used for its side effect of generating relevant
##' parts of the XML file.
##' @author François Michonneau
randomizeTraits <- function(traitFile, output) {

    tf <-  read.table(file=traitFile, header=T, stringsAsFactors=FALSE)
    tf[, 2] <- sample(tf[, 2])

    outVec1 <- character(nrow(tf))
    for (i in 1:nrow(tf)) {
        outVec1[i] <- paste("\t\t<taxon id=\"", tf[i,1], "\">\n\t\t\t<attr name=\"species\">\n\t\t\t\t",
                           tf[i,2], "\n\t\t\t</attr>\n\t\t</taxon>", sep="")
    }


    allSpp <- unique(tf[, 2])
    outVec2 <- character(2 * length(allSpp) + nrow(tf))
    k <- 1
    for (i in 1:length(allSpp)) {
        tmpDt <- subset(tf, species == allSpp[i])
        outVec2[k] <- paste("\t\t<sp id=\"", allSpp[i], "\">", sep="")
        k <- k + 1
        for (j in 1:nrow(tmpDt)) {
          outVec2[k] <- paste ("\t\t\t<taxon idref=\"", tmpDt[j, 1], "\"/>", sep="")
          k <- k + 1
        }
        outVec2[k] <- "\t\t</sp>"
        k <- k + 1
    }
    
    cat(outVec1, "\n\n----------------\n\n", outVec2, file="/tmp/xx.txt", sep="\n")
}

