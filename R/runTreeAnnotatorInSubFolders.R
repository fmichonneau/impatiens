
##' Given a folder, run treeannotator in all subfolders
##'
##' .. content for \details{} ..
##' @title 
##' @param path 
##' @return 
##' @author Francois Michonneau
runTreeAnnotatorInSubFolder <- function(path,
                                        treeAnnotatorPath="~/Software/BEASTv1.8.0/bin/./treeannotator",
                                        treeAnnotatorOpts="-heights ca -burnin 1000",
                                        ...) {

    inFiles <- list.files(pattern="trees$", path=path, full.names=TRUE, recursive=TRUE)
    outFiles <- gsub("trees$", "tre.nex", inFiles)
    cmds <- paste(treeAnnotatorPath, treeAnnotatorOpts, inFiles, outFiles)
    for (i in 1:length(inFiles)) {
        system(cmds[i])
    }    
}
