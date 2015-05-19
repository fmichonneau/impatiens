check_file <- function(file, verbose) {
    if (!file.exists(file)) {
        stop(file, " wasn't created.")
    } else {
        if (verbose) message(file, " succesfully created.")
    }
}

create_dir <- function(file, verbose) {
    if (!file.exists(dirname(file))) {
        if (verbose) {
            message("Create: ", dirname(file))
        }
        dir.create(dirname(file), recursive = TRUE)
    }
}

make_csv <- function(obj, file, ...,  verbose = TRUE) {
    on.exit(check_file(file, verbose = verbose))
    create_dir(file = file, verbose = verbose)
    if (verbose) {
        message("Creating csv file: ", file)
    }
    write.csv(obj, file = file, row.names = FALSE, ...)
}
