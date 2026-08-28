check_log_transformation <- function(expression_data) {
    max_val <- max(expression_data, na.rm = TRUE)
    min_val <- min(expression_data, na.rm = TRUE)
    cat("Expression data summary:\n")
    cat("  Min value:", min_val, "\n")
    cat("  Max value:", max_val, "\n")
    cat("  Range:", max_val - min_val, "\n")
    
    # Check if data appears to be already log2-transformed
    if (max_val < 50 && min_val > 0) {
        cat("  Assessment: Data appears to be already log2-transformed\n")
        cat("  Recommendation: No additional log2 transformation", 
            "needed\n")
        needs_log <- FALSE
    } else if (max_val > 100) {
        cat("  Assessment: Data appears to be in raw scale", 
            "(not log-transformed)\n")
        cat("  Recommendation: Apply log2 transformation\n")
        needs_log <- TRUE
    } else {
        cat("  Assessment: Unclear - manual inspection recommended\n")
        cat("  Recommendation: Check data distribution and source", 
            "documentation\n")
        needs_log <- NA
    }
    
    # Apply transformation based on assessment
    if (isTRUE(needs_log)) {
        cat("\nApplying log2 transformation...\n")
        # Add small constant to avoid log(0)
        expression_data <- log2(expression_data + 1)
        cat("Transformation complete.\n")
    } else if (isFALSE(needs_log)) {
        cat("\nNo transformation needed - using data as is.\n")
    } else {
        cat("\nPlease manually verify if transformation is needed.\n")
    }
    
    return(expression_data)
}
