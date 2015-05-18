choose_dev <- function() {
    if (identical(Sys.getenv("DEV_TYPE"), "tikz")) {
        return("tikz")
    } else {
        return("pdf")
    }
}
