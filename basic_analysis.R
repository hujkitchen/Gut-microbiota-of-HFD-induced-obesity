#=============================================================================
# BASIC LONGITUDINAL MICROBIOME ANALYSIS PIPELINE
# Direct data extraction from TreeSummarizedExperiment
#=============================================================================

# ============================================================================
# SETUP AND LIBRARIES
# ============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("ggplot2", "patchwork", "dplyr", "tidyr", "tibble", 
                       "lme4", "lmerTest", "pheatmap", "RColorBrewer", 
                       "ggpubr", "viridis", "ggrepel", "vegan", "broom", "purrr")

for(pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

if (!require("maaslin3", quietly = TRUE)) {
    BiocManager::install("maaslin3")
    library(maaslin3)
}

# Set options
options(stringsAsFactors = FALSE)
set.seed(123)
theme_set(theme_bw(base_size = 12))

# Colors
diet_colors <- c("ND" = "#2E86AB", "HFD" = "#A23B72")
phase_colors <- c("Acute" = "#440154", "Early" = "#3B528B", 
                  "Mid" = "#21908C", "Late" = "#5DC863")

# Create output directory
output_dir <- "microbiome_analysis_results"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

# ============================================================================
# DIRECT DATA EXTRACTION FROM TSE
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("EXTRACTING DATA FROM TreeSummarizedExperiment\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Extract count matrix directly
cat("Extracting count matrix...\n")
counts_raw <- assay(tse, "counts")

# Check what we have
cat("Raw counts class:", class(counts_raw), "\n")
cat("Raw counts dimensions:", dim(counts_raw), "\n")

# Convert to plain matrix if it's a special class
if (inherits(counts_raw, "dgCMatrix")) {
    cat("Converting from sparse matrix to dense matrix...\n")
    counts_raw <- as.matrix(counts_raw)
}

# Ensure numeric
count_matrix <- matrix(as.numeric(counts_raw), 
                       nrow = nrow(counts_raw), 
                       ncol = ncol(counts_raw),
                       dimnames = dimnames(counts_raw))

cat("Count matrix class after conversion:", class(count_matrix), "\n")
cat("Count matrix type:", typeof(count_matrix), "\n")

# Check for NA/NaN
cat("NA values:", sum(is.na(count_matrix)), "\n")
cat("NaN values:", sum(is.nan(count_matrix)), "\n")
cat("Infinite values:", sum(is.infinite(count_matrix)), "\n")

# Replace any NA/NaN with 0
count_matrix[is.na(count_matrix) | is.nan(count_matrix)] <- 0
count_matrix[is.infinite(count_matrix)] <- 0

# Extract metadata
cat("\nExtracting sample metadata...\n")
metadata <- as.data.frame(colData(tse))

cat("Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

# Check column types
cat("\nMetadata column types:\n")
print(sapply(metadata, class))

# Extract taxonomic information
cat("\nExtracting taxonomic information...\n")
taxonomy <- as.data.frame(rowData(tse))
cat("Taxonomy columns:", paste(colnames(taxonomy), collapse = ", "), "\n")
cat("Taxonomy dimensions:", dim(taxonomy), "\n")

# ============================================================================
# PREPARE METADATA - FIXED
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PREPARING METADATA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Check weight column
cat("Weight column type:", class(metadata$weight), "\n")
cat("First few weight values:\n")
print(head(metadata$weight))

# Clean and prepare metadata - ensure all numeric columns are properly converted
metadata_clean <- metadata

# Convert timepoint and extract day
metadata_clean$Day <- as.numeric(gsub(".*D", "", as.character(metadata_clean$timepiont)))

# Create timepoint factor
timepoint_levels <- paste0("D", sort(unique(metadata_clean$Day)))
metadata_clean$timepoint <- factor(metadata_clean$timepiont, 
                                    levels = timepoint_levels,
                                    ordered = TRUE)

# Create phases
metadata_clean$Phase <- cut(metadata_clean$Day,
                             breaks = c(0, 7, 21, 42, 71),
                             labels = c("Acute", "Early", "Mid", "Late"))

# Convert weight to numeric if needed
if (!is.numeric(metadata_clean$weight)) {
    cat("Converting weight to numeric...\n")
    metadata_clean$weight <- as.numeric(as.character(metadata_clean$weight))
}

# Check weight after conversion
cat("Weight after conversion - type:", class(metadata_clean$weight), "\n")
cat("Weight summary:\n")
print(summary(metadata_clean$weight))

# Create scaled variables safely
metadata_clean$weight_scaled <- as.numeric(scale(metadata_clean$weight))
metadata_clean$day_scaled <- as.numeric(scale(metadata_clean$Day))

# Ensure factors
metadata_clean$group1 <- factor(metadata_clean$group1)
metadata_clean$cage <- factor(metadata_clean$cage)

# Add SampleID as rownames if not present
if (!"SampleID" %in% colnames(metadata_clean)) {
    metadata_clean$SampleID <- colnames(count_matrix)
}

cat("\nSample distribution by diet and timepoint:\n")
print(table(metadata_clean$group1, metadata_clean$timepoint))

cat("\nSample distribution by diet and phase:\n")
print(table(metadata_clean$group1, metadata_clean$Phase))

cat("\nCage distribution:\n")
print(table(metadata_clean$cage, metadata_clean$group1))

# ============================================================================
# SEQUENCING DEPTH ANALYSIS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SEQUENCING DEPTH ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate depth using apply
depth_values <- apply(count_matrix, 2, sum, na.rm = TRUE)

depth_data <- data.frame(
    Sample = colnames(count_matrix),
    Depth = depth_values,
    Day = metadata_clean$Day,
    Diet = metadata_clean$group1,
    Cage = metadata_clean$cage,
    Phase = metadata_clean$Phase
)

cat("Sequencing depth summary:\n")
cat("Mean depth:", round(mean(depth_values)), "\n")
cat("Min depth:", min(depth_values), "\n")
cat("Max depth:", max(depth_values), "\n")

# Depth plots
p_depth1 <- ggplot(depth_data, aes(x = Day, y = Depth, color = Diet)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(aes(fill = Diet), method = "loess", se = TRUE, alpha = 0.2) +
    labs(x = "Day", y = "Sequencing Depth", 
         title = "Sequencing Depth Over Time") +
    scale_color_manual(values = diet_colors) +
    scale_fill_manual(values = diet_colors) +
    theme(legend.position = "bottom")

p_depth2 <- ggplot(depth_data, aes(x = Phase, y = Depth, fill = Diet)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    labs(x = "Phase", y = "Sequencing Depth", 
         title = "Sequencing Depth by Phase") +
    scale_fill_manual(values = diet_colors) +
    scale_y_log10()

depth_combined <- p_depth1 | p_depth2
ggsave("sequencing_depth_analysis.pdf", depth_combined, width = 12, height = 6)
print(depth_combined)

# ============================================================================
# FILTERING
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FILTERING LOW PREVALENCE TAXA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Prevalence filtering
prev_cutoff <- 0.1
min_samples <- ceiling(ncol(count_matrix) * prev_cutoff)
keep_taxa <- apply(count_matrix > 0, 1, sum) >= min_samples

cat("Prevalence cutoff:", prev_cutoff, "\n")
cat("Minimum samples required:", min_samples, "\n")
cat("Taxa before filtering:", nrow(count_matrix), "\n")
cat("Taxa after filtering:", sum(keep_taxa), "\n")

# Filter count matrix
count_matrix_filtered <- count_matrix[keep_taxa, ]
taxonomy_filtered <- taxonomy[keep_taxa, ]

# ============================================================================
# NORMALIZATION
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("NORMALIZATION\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate relative abundances
col_sums <- apply(count_matrix_filtered, 2, sum)
rel_abundance <- sweep(count_matrix_filtered, 2, col_sums, "/")

cat("Relative abundance column sums (should be 1):\n")
print(summary(colSums(rel_abundance)))

# ============================================================================
# AGGREGATE TO GENUS LEVEL
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("AGGREGATING TO GENUS LEVEL\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Check if genus column exists
if ("genus" %in% colnames(taxonomy_filtered)) {
    cat("Genus column found. Aggregating...\n")
    
    # Get genus for each row
    genus_labels <- taxonomy_filtered$genus
    genus_labels[is.na(genus_labels)] <- "Unknown"
    
    # Aggregate by genus
    unique_genera <- unique(genus_labels)
    cat("Number of unique genera:", length(unique_genera), "\n")
    
    # Create genus-level count matrix
    genus_counts <- matrix(0, nrow = length(unique_genera), ncol = ncol(count_matrix_filtered))
    rownames(genus_counts) <- unique_genera
    colnames(genus_counts) <- colnames(count_matrix_filtered)
    
    # Use a loop or aggregate
    for (g in unique_genera) {
        idx <- which(genus_labels == g)
        if (length(idx) > 1) {
            genus_counts[g, ] <- colSums(count_matrix_filtered[idx, , drop = FALSE])
        } else {
            genus_counts[g, ] <- count_matrix_filtered[idx, ]
        }
    }
    
    cat("Genus count matrix dimensions:", dim(genus_counts), "\n")
    
    # Genus-level relative abundance
    genus_rel <- sweep(genus_counts, 2, colSums(genus_counts), "/")
    
} else {
    cat("Genus column not found. Using ASV level.\n")
    genus_counts <- count_matrix_filtered
    genus_rel <- rel_abundance
    rownames(genus_counts) <- paste0("ASV", 1:nrow(genus_counts))
}

# ============================================================================
# ALPHA DIVERSITY
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ALPHA DIVERSITY ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate alpha diversity
count_matrix_t <- t(count_matrix_filtered)

alpha_df <- metadata_clean
alpha_df$shannon <- diversity(count_matrix_t, index = "shannon")
alpha_df$simpson <- diversity(count_matrix_t, index = "simpson")
alpha_df$observed <- specnumber(count_matrix_t)

cat("Alpha diversity summary by diet:\n")
alpha_summary <- alpha_df %>%
    group_by(group1) %>%
    summarise(
        Mean_Shannon = mean(shannon, na.rm = TRUE),
        SD_Shannon = sd(shannon, na.rm = TRUE),
        Mean_Observed = mean(observed, na.rm = TRUE),
        SD_Observed = sd(observed, na.rm = TRUE),
        .groups = 'drop'
    )
print(alpha_summary)

# Alpha diversity plots
p_alpha1 <- ggplot(alpha_df, aes(x = Day, y = shannon, color = group1)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(aes(fill = group1), method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
    labs(x = "Day", y = "Shannon Diversity", 
         title = "Shannon Diversity Trajectory") +
    scale_color_manual(values = diet_colors, name = "Diet") +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    theme(legend.position = "bottom")

p_alpha2 <- ggplot(alpha_df, aes(x = Phase, y = shannon, fill = group1)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    labs(x = "Phase", y = "Shannon Diversity", 
         title = "Shannon Diversity by Phase") +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    theme(legend.position = "none")

p_alpha3 <- ggplot(alpha_df, aes(x = Day, y = observed, color = group1)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(aes(fill = group1), method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
    labs(x = "Day", y = "Observed Richness", 
         title = "Species Richness Trajectory") +
    scale_color_manual(values = diet_colors, name = "Diet") +
    scale_fill_manual(values = diet_colors, name = "Diet") +
    theme(legend.position = "none")

alpha_combined <- (p_alpha1 | p_alpha3) / p_alpha2 +
    plot_annotation(tag_levels = 'A')

ggsave("alpha_diversity_plots.pdf", alpha_combined, width = 12, height = 10)
print(alpha_combined)

# Statistical test
cat("\nAlpha Diversity Statistical Test:\n")
tryCatch({
    shannon_lm <- lm(shannon ~ group1 * Day + weight_scaled + cage, data = alpha_df)
    cat("\nShannon Diversity - ANOVA:\n")
    print(anova(shannon_lm))
}, error = function(e) {
    cat("Linear model failed:", e$message, "\n")
    # Try simpler model
    shannon_lm <- lm(shannon ~ group1 + Day + cage, data = alpha_df)
    cat("\nSimpler model ANOVA:\n")
    print(anova(shannon_lm))
})

# ============================================================================
# BETA DIVERSITY
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("BETA DIVERSITY ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate Bray-Curtis distance
bray_dist <- vegdist(t(rel_abundance), method = "bray")

# PCoA
pcoa_res <- cmdscale(bray_dist, k = 5, eig = TRUE)
var_exp <- round(pcoa_res$eig[1:3] / sum(pcoa_res$eig) * 100, 1)

pcoa_df <- data.frame(
    PC1 = pcoa_res$points[,1],
    PC2 = pcoa_res$points[,2],
    PC3 = pcoa_res$points[,3],
    alpha_df
)

# PCoA plots
p_beta1 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group1)) +
    geom_point(aes(size = Day), alpha = 0.7) +
    stat_ellipse(aes(group = group1), level = 0.95, alpha = 0.3) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"),
         title = "Bray-Curtis PCoA: Diet Effect") +
    scale_color_manual(values = diet_colors, name = "Diet") +
    scale_size_continuous(name = "Day")

p_beta2 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Day)) +
    geom_point(size = 2.5, alpha = 0.7) +
    facet_wrap(~ group1) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"),
         title = "Bray-Curtis PCoA: Temporal Trajectory") +
    scale_color_viridis_c(option = "plasma", name = "Day")

beta_combined <- p_beta1 / p_beta2 +
    plot_annotation(tag_levels = 'A')

ggsave("beta_diversity_pcoa.pdf", beta_combined, width = 10, height = 10)
print(beta_combined)

# PERMANOVA
cat("\nPERMANOVA Results:\n")
tryCatch({
    permanova_res <- adonis2(bray_dist ~ group1 * Phase + weight_scaled, 
                              data = alpha_df, permutations = 999)
    print(permanova_res)
}, error = function(e) {
    cat("PERMANOVA failed:", e$message, "\n")
    permanova_res <- adonis2(bray_dist ~ group1 * Phase, 
                              data = alpha_df, permutations = 999)
    print(permanova_res)
})

# ============================================================================
# DIFFERENTIAL ABUNDANCE WITH MaAsLin3
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("DIFFERENTIAL ABUNDANCE ANALYSIS WITH MaAsLin3\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Prepare metadata for MaAsLin3
metadata_maaslin <- data.frame(
    group1 = factor(metadata_clean$group1),
    Day = metadata_clean$Day,
    Phase = factor(metadata_clean$Phase),
    cage = factor(metadata_clean$cage),
    weight_scaled = metadata_clean$weight_scaled,
    day_scaled = metadata_clean$day_scaled,
    row.names = metadata_clean$SampleID
)

# Use genus-level counts
genus_counts_df <- as.data.frame(genus_counts)

cat("Running MaAsLin3 for diet effect...\n")
cat("Input dimensions:", dim(genus_counts_df), "\n")

tryCatch({
    maaslin_results <- maaslin3(
        input_data = genus_counts_df,
        input_metadata = metadata_maaslin,
        output = "maaslin3_results",
        min_abundance = 0.0001,
        min_prevalence = 0.1,
        normalization = "TSS",
        transform = "LOG",
        analysis_method = "LM",
        fixed_effects = c("group1", "day_scaled", "weight_scaled"),
        random_effects = c("cage"),
        reference = c("group1,ND"),
        max_significance = 0.05,
        plot_heatmap = TRUE,
        plot_scatter = TRUE,
        cores = 4
    )
    
    # Extract significant results
    sig_results <- maaslin_results$results %>%
        filter(metadata == "group1", qval < 0.05) %>%
        arrange(qval)
    
    cat("\nSignificant diet-associated features:", nrow(sig_results), "\n")
    
    if (nrow(sig_results) > 0) {
        write.csv(sig_results, "maaslin3_diet_effects.csv", row.names = FALSE)
        cat("\nTop 20 diet-associated features:\n")
        print(head(sig_results[, c("feature", "coef", "pval", "qval")], 20))
    }
}, error = function(e) {
    cat("MaAsLin3 with random effects failed:", e$message, "\n")
    cat("\nTrying without random effects...\n")
    
    tryCatch({
        maaslin_results <- maaslin3(
            input_data = genus_counts_df,
            input_metadata = metadata_maaslin,
            output = "maaslin3_results_noRE",
            min_abundance = 0.0001,
            min_prevalence = 0.1,
            normalization = "TSS",
            transform = "LOG",
            analysis_method = "LM",
            fixed_effects = c("group1", "day_scaled", "weight_scaled"),
            reference = c("group1,ND"),
            max_significance = 0.05,
            cores = 4
        )
        
        sig_results <- maaslin_results$results %>%
            filter(metadata == "group1", qval < 0.05) %>%
            arrange(qval)
        
        cat("\nSignificant diet-associated features (no RE):", nrow(sig_results), "\n")
        
        if (nrow(sig_results) > 0) {
            write.csv(sig_results, "maaslin3_diet_effects_noRE.csv", row.names = FALSE)
        }
    }, error = function(e2) {
        cat("MaAsLin3 without random effects also failed:", e2$message, "\n")
    })
})

# ============================================================================
# SAVE RESULTS AND SUMMARY
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Save workspace
save.image("complete_analysis_workspace.RData")

# Save filtered data
saveRDS(list(
    count_matrix = count_matrix_filtered,
    metadata = metadata_clean,
    taxonomy = taxonomy_filtered,
    genus_counts = genus_counts,
    genus_rel = genus_rel,
    alpha_diversity = alpha_df
), file = "processed_data.rds")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("All results saved in:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - sequencing_depth_analysis.pdf\n")
cat("  - alpha_diversity_plots.pdf\n")
cat("  - beta_diversity_pcoa.pdf\n")
cat("  - processed_data.rds\n")
cat("  - complete_analysis_workspace.RData\n")

# Print summary
cat("\nSummary:\n")
cat("  Total samples:", ncol(count_matrix), "\n")
cat("  Filtered taxa:", nrow(count_matrix_filtered), "\n")
cat("  Genera:", nrow(genus_counts), "\n")
cat("  Timepoints:", paste(levels(metadata_clean$timepoint), collapse = ", "), "\n")