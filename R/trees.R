load_impTree <- function(filename) {
    impTree <- phyloch::read.beast(file="data/allimpatiens_strict.tree.nex")
    impTree
}

load_impTree4 <- function(tree) {
    return(as(tree, "phylo4"))
}


drop_tip_from_raxml_trees <- function(tree_file, outfile, to.drop = "S0213") {
    allRaxmlBS <- ape::read.tree(file=tree_file)
    noOutgroup <- lapply(allRaxmlBS, function(tr) drop.tip(tr, to.drop))
    class(noOutgroup) <- "multiPhylo"
    ape::write.tree(noOutgroup, file = outfile)
}

annotate_beast_with_bootstrap <- function(beast_tree, raxml_tree,
                                          raxml_path = "~/Software/RAxML-8.0.1/./raxmlHPC-PTHREADS-SSE3") {
    rxmlCmd <- paste(raxml_path, "-m GTRGAMMA",
                     "-p 12345 -f b -t", beast_tree,
                     "-z", raxml_tree, "-T8 -n annotateBEASTtree")
    system(rxmlCmd)
    system("mv *.annotateBEASTtree data/")
    stopifnot(identical(list.files(pattern="annotateBEASTtree$", path="data/"),
                        c("RAxML_bipartitions.annotateBEASTtree",
                          "RAxML_bipartitionsBranchLabels.annotateBEASTtree",
                          "RAxML_info.annotateBEASTtree")))
}

get_prop_clades <- function(raxml_trees, rds_file){
    noOutgroup <- ape::read.tree(file=raxml_trees)
    ppNoOutGroup <- ape::prop.part(noOutgroup)
    pcNoOutGroup <- ape::prop.clades(impTree, part=ppNoOutGroup)
    saveRDS(pcNoOutGroup, file=rds_file)
}

load_per_locus_raxml_tree_files <- function() {
    loci <- c("c0036", "c0775", "H3a", "mtDNA", "rDNA")
    raxFiles <- file.path("data/raxml_perlocus", paste0("RAxML_bipartitions.", loci))
    stopifnot(all(sapply(raxFiles, file.exists)))
    raxFiles
}
