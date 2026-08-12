#' Read .RData network files and Combine them
#'
#' This function reads individual .RData network files from specified file 
#' paths and combines them together if there is more than one .RData file. 
#'
#' @param network_files vector of file paths to the network .RData files
#' @return A combined network object containing all the networks from the
#'         specified .RData files.
#' @export
#'
#' 
read_networks <- function(network_files) {
    # Ensure we have file paths
    if (length(network_files) == 0) {
        stop("No network files provided")
    }
    # Check if files exist
    missing_files <- network_files[!file.exists(network_files)]
    if (length(missing_files) > 0) {
        stop(sprintf("The following files do not exist:\n%s", 
                    paste(missing_files, collapse = "\n")))
    }
    cat("Following files will be processed:", "\n")
    print(basename(network_files))
    networks_all <- NULL
    for (i in seq_along(network_files)) {
        cat("Loading file", basename(network_files[i]), "\n")
        load(network_files[i], data <- new.env())
        networks <- data[["net"]]
        networks_all <- cbind(networks_all, networks)
        rm(networks); gc()
    }
    return(networks_all)
}
#' Quantile Normalize Networks
#'
#' This function reads network files from specified file paths and applies
#' quantile normalization to the combined network data.
#'
#' @param network_files vector of file paths to the network .RData files
#' @return A matrix representing the quantile normalized network.
#' @export
#'
quantile_normalize_net <- function(network_files) {
    cat("Reading in networks", "\n")
    net <- read_networks(network_files)
    cat("Normalizing network", "\n")
    net_norm <- normalize.quantiles(as.matrix(net))
    colnames(net_norm) <- colnames(net)
    return(net_norm)
    rm(net)
    gc()
}