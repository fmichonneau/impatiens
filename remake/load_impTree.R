load_impTree <- function(filename) {
    impTree <- read.beast(file="data/allimpatiens_strict.tree.nex")
    impTree
}

load_impTree4 <- function(tree) {
    return(as(tree, "phylo4"))
}
