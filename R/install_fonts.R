install_fonts <- function() {
    dir.create("~/.texmf-var/luatex-cache/generic/fonts/otf")
    system("git clone git@github.com:fmichonneau/latex-fonts.git ~/.texmf-var/luatex-cache/generic/fonts/otf")
    stopifnot(length(list.files(path = "~/.texmf-var/luatex-cache/generic/fonts/otf", pattern = "lua$|luc$")) > 1)
}
