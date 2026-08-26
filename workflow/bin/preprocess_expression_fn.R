#' Create PCA Plot with Outlier Detection

create_pca_outlier_plot <- function(
    data_matrix,
    sample_names = NULL,
    covariate = NULL,
    outlier_method = "distance",
    outlier_threshold = 2,
    scale_data = TRUE,
    center_data = TRUE,
    plot_title = "PCA Plot",
    color_palette = c("FALSE" = "black", "TRUE" = "red"),
    point_size = 3,
    point_alpha = 0.7,
    label_outliers = TRUE) {

    sample_names <- colnames(data_matrix)

    # Perform PCA (transpose so samples are rows)
    cat("=== PCA ANALYSIS ===\n")
    pca_res <- prcomp(
        t(data_matrix),
        scale. = scale_data,
        center = center_data
    )

    # Create PCA data frame
    pca_data <- data.frame(
        Sample = sample_names,
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2],
        stringsAsFactors = FALSE
    )

    # Add covariate if provided
    if (!is.null(covariate)) {
        if (length(covariate) != ncol(data_matrix)) {
            stop(
                "Length of covariate vector must match number of samples ",
                "in data_matrix"
            )
        }
        pca_data$Covariate <- covariate
    }

    # Calculate variance explained
    var_explained <- summary(pca_res)$importance[2, ]
    pc1_var <- round(var_explained[1] * 100, 1)
    pc2_var <- round(var_explained[2] * 100, 1)

    pca_center <- c(
        mean(pca_data$PC1),
        mean(pca_data$PC2)
    )

    pca_data$distance_from_center <- sqrt(
        (pca_data$PC1 - pca_center[1])^2 +
            (pca_data$PC2 - pca_center[2])^2
    )

    # Define outliers
    dist_threshold <- mean(pca_data$distance_from_center) +
        outlier_threshold * sd(pca_data$distance_from_center)

    pca_data$is_outlier <-
        pca_data$distance_from_center > dist_threshold

    # Report outlier information
    n_outliers <- sum(pca_data$is_outlier)
    outlier_samples <- pca_data[pca_data$is_outlier, ]

    # Create PCA plot
    cat("\n=== CREATING PLOT ===\n")

    if (!is.null(covariate)) {
        # Color by covariate
        p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
            geom_point(
                aes(color = Covariate),
                size = point_size,
                alpha = point_alpha
            ) +
            theme_minimal() +
            labs(
                title = plot_title,
                x = paste0("PC1 (", pc1_var, "%)"),
                y = paste0("PC2 (", pc2_var, "%)"),
                color = "Covariate"
            )
    } else {
        # Color by outlier status (default behavior)
        p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
            geom_point(
                aes(color = is_outlier),
                size = point_size,
                alpha = point_alpha
            ) +
            scale_color_manual(
                values = color_palette,
                labels = c("Normal", "Outlier"),
                name = "Sample Type"
            ) +
            theme_minimal() +
            labs(
                title = plot_title,
                x = paste0("PC1 (", pc1_var, "%)"),
                y = paste0("PC2 (", pc2_var, "%)")
            )
    }

    # Add outlier labels if requested
    # (only when coloring by outliers)
    if (is.null(covariate) && label_outliers && n_outliers > 0) {
        p <- p +
            geom_text_repel(
                data = pca_data[pca_data$is_outlier, ],
                aes(label = Sample),
                size = 3,
                box.padding = 0.5,
                point.padding = 0.3,
                segment.color = color_palette["TRUE"],
                segment.alpha = 0.6
            )
    }

    print(p)
    print(outlier_samples)
}