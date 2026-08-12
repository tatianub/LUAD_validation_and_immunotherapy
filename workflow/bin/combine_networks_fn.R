#' Read and combine LIONESS network files
#'
#' Reads individual .txt network files within a specified range and combines
#' them into a single data frame with samples as columns.
#'
#' @param networks_cancer Data frame with \code{file_path} and \code{tcga_id}
#'   columns containing network file paths and sample identifiers
#' @param start Integer. First network index to read
#' @param end Integer. Last network index to read
#'
#' @return Data frame with combined networks where rows are TF-gene edges
#'   and columns are samples (named by TCGA IDs)
#'
#' @export

read_networks_txt <- function(networks_cancer, start, end) {
    file_paths <- networks_cancer$file_path[start:end]
    cat("Reading", length(file_paths), "network files...\n")
    
    # Check if files exist
    missing_files <- file_paths[!file.exists(file_paths)]
    if (length(missing_files) > 0) {
        cat("ERROR: Missing files:\n")
        cat(paste(missing_files[1:min(5, length(missing_files))], collapse = "\n"), "\n")
        if (length(missing_files) > 5) cat("... and", length(missing_files) - 5, "more\n")
        stop("Cannot find network files")
    }
    
    datalist <- 
        lapply(file_paths, function(x) fread(x))
    datalist <- lapply(datalist, function(x) melt(x))
    datalist <- map(datalist, ~ (.x %>% select(-variable)))
    net <- dplyr::bind_cols(datalist)
    colnames(net) <- networks_cancer$sample_id[start:end]
    return(net)
}