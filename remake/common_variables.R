load_esuList <- function() {
    c("ESU1", "ESU2", "ESU3", "gracilis", "tiger", "tigerRedSea", "Medit", "WA",
      "Gala", "EP", "Hawaii", "Wpac", "RedSea")
}

load_gmycFactors <- function() {
    gmycFactors <- expand.grid(c("strict", "relaxed"),
                               c("yule", "coalexp", "coalcst"))
    gmycFactors <- apply(gmycFactors, 1, paste0, collapse="_")
    gmycFactors <- c(paste0("20140422.allSeq_", gmycFactors),
                     paste0("20140514.noDup_", gmycFactors))
    gmycFactors
}


load_impPal <- function() {
    impPal <- c(
        "#89B151", # tenenbaum bright green = Medit
        "#AF420A", # life aquatic dark orange = ESU3_Deep
        "#FD6467", # grand budapest hotel pink = ESU3_PNG
        "#01abe9", # life aquatic bright blue = ESU1
        "#1b346c", # life aquatic dark blue = ESU3
        "#f54b1a", # life aquatic bright orange = gracilis
        "#e5c39e", # life aquatic tan = Hawaii
        "#c3ced0", # life aquatic light blue = RedSea
        "#EBCC2A", # life aquatic yellow = Wpac
        "#446455", # tenenbaum dark green = Gala
        "#4CBDD5", # darjeeling light blue = WA
        "#FF0000", # darjeeling red = tigerRedSea
        "#F8B0BB", # pink = tiger
        "#B6F9B1", # another green = ESU1_Lizard
        "#FFF196", # mr fox yellow = EP
        "#8061AD" # rushmore purple = ESU2
        ##"#020202"  # not quite black
    )
    names(impPal) <- c("Medit", "ESU3_Deep", "ESU3_PNG", "ESU1", "ESU3",
                       "gracilis", "Hawaii", "RedSea", "Wpac", "Gala", "WA",
                       "tigerRedSea", "tiger", "ESU1_Lizard", "EP", "ESU2")
    impPal
}

get_nuclear_loci <- function() {
    c("c0036", "c0775", "H3a", "ITS", "LSU")
}

get_mt_loci <- function() {
    c("16S", "ATP6", "COI")
}

get_all_loci <- function() {
    c(get_nuclear_loci(), get_mt_loci())
}
