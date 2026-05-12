#=============================================================================
# ANCOM-BC2 DIFFERENTIAL ABUNDANCE ANALYSIS
# Using correct function arguments
#=============================================================================

# ============================================================================
# INSTALL AND LOAD ANCOM-BC2
# ============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("ANCOMBC", quietly = TRUE)) {
    BiocManager::install("ANCOMBC")
}

library(ANCOMBC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# Check ANCOM-BC2 function arguments
cat("\nANCOM-BC2 function arguments:\n")
print(args(ancombc2))

# ============================================================================
# PREPARE DATA - Convert TSE to phyloseq object
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PREPARING DATA FOR ANCOM-BC2\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ANCOM-BC2 works with phyloseq or TreeSummarizedExperiment objects
# Let's create a proper TreeSummarizedExperiment object

# Prepare count matrix (ensure it's integer)
genus_counts_int <- round(genus_counts)
mode(genus_counts_int) <- "integer"

# Prepare metadata
sample_meta <- data.frame(
    group1 = factor(metadata_clean$group1),
    Day = metadata_clean$Day,
    Phase = factor(metadata_clean$Phase),
    cage = factor(metadata_clean$cage),
    weight_scaled = metadata_clean$weight_scaled,
    day_scaled = metadata_clean$day_scaled,
    row.names = colnames(genus_counts)
)

# Create a simple TreeSummarizedExperiment
library(SummarizedExperiment)

tse_genus <- TreeSummarizedExperiment(
    assays = list(counts = genus_counts_int),
    colData = sample_meta,
    rowData = data.frame(genus = rownames(genus_counts_int),
                         row.names = rownames(genus_counts_int))
)

cat("TSE object created:\n")
cat("Dimensions:", dim(tse_genus), "\n")

# ============================================================================
# ANCOM-BC2 ANALYSIS 1: DIET EFFECT
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANCOM-BC2: DIET EFFECT (HFD vs ND)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Correct ANCOM-BC2 syntax
# The function uses 'data' parameter which is a TSE or phyloseq object
# and 'fix_formula' for the model formula

ancom_diet <- tryCatch({
    ancombc2(
        data = tse_genus,                    # TSE object
        assay_name = "counts",               # Use assay_name instead of assay.type
        tax_level = "genus",                 # Taxonomic level
        fix_formula = "group1 + day_scaled + weight_scaled",
        rand_formula = "(1 | cage)",
        p_adj_method = "holm",
        prv_cut = 0.1,
        lib_cut = 1000,
        group = "group1",
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        n_cl = 4
    )
}, error = function(e) {
    cat("ANCOM-BC2 with random effects failed:", e$message, "\n")
    cat("Trying without random effects...\n\n")
    
    tryCatch({
        ancombc2(
            data = tse_genus,
            assay_name = "counts",
            tax_level = "genus",
            fix_formula = "group1 + day_scaled + weight_scaled",
            p_adj_method = "holm",
            prv_cut = 0.1,
            lib_cut = 1000,
            group = "group1",
            struc_zero = TRUE,
            neg_lb = TRUE,
            alpha = 0.05,
            n_cl = 4
        )
    }, error = function(e2) {
        cat("ANCOM-BC2 still failed:", e2$message, "\n\n")
        cat("Trying with minimal parameters...\n")
        
        ancombc2(
            data = tse_genus,
            assay_name = "counts",
            tax_level = "genus",
            fix_formula = "group1",
            p_adj_method = "holm",
            prv_cut = 0.1,
            lib_cut = 1000,
            group = "group1",
            struc_zero = TRUE,
            neg_lb = TRUE
        )
    })
})

cat("\nANCOM-BC2 analysis completed successfully!\n")

# ============================================================================
# EXTRACT AND PROCESS RESULTS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PROCESSING ANCOM-BC2 RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Check the structure of the results
cat("Result object names:\n")
print(names(ancom_diet))

# Extract results from the output
res <- ancom_diet$res

cat("\nResult columns:\n")
print(colnames(res))

# Display first few rows
cat("\nFirst 5 rows of results:\n")
print(head(res, 5))

# ============================================================================
# IDENTIFY CORRECT COLUMN NAMES
# ============================================================================

# Find group1-related columns
group1_cols <- grep("group1", colnames(res), value = TRUE)
cat("\nGroup1 related columns:", paste(group1_cols, collapse = ", "), "\n")

# Extract diet effect results
lfc_col <- grep("^lfc_group1", colnames(res), value = TRUE)
se_col <- grep("^se_group1", colnames(res), value = TRUE)
p_col <- grep("^p_group1", colnames(res), value = TRUE)
q_col <- grep("^q_group1", colnames(res), value = TRUE)
diff_col <- grep("^diff_group1", colnames(res), value = TRUE)

cat("\nIdentified columns:\n")
cat("  LFC:", lfc_col, "\n")
cat("  SE:", se_col, "\n")
cat("  P-value:", p_col, "\n")
cat("  Q-value:", q_col, "\n")
cat("  Diff:", diff_col, "\n")

# ============================================================================
# CREATE RESULTS DATAFRAME
# ============================================================================

ancom_results <- data.frame(
    taxon = res$taxon,
    stringsAsFactors = FALSE
)

# Add columns if they exist
if(length(lfc_col) > 0) ancom_results$lfc <- res[[lfc_col]]
if(length(se_col) > 0) ancom_results$se <- res[[se_col]]
if(length(p_col) > 0) ancom_results$pval <- res[[p_col]]
if(length(q_col) > 0) ancom_results$qval <- res[[q_col]]
if(length(diff_col) > 0) ancom_results$diff <- res[[diff_col]]

# Add W statistic if present
w_col <- grep("^W_group1", colnames(res), value = TRUE)
if(length(w_col) > 0) ancom_results$W <- res[[w_col]]

# Summary statistics
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("DIET EFFECT RESULTS SUMMARY\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")
cat("Total genera tested:", nrow(ancom_results), "\n")

if("diff" %in% colnames(ancom_results)) {
    cat("Differentially abundant genera:", sum(ancom_results$diff, na.rm = TRUE), "\n")
    
    sig_ancom <- ancom_results %>%
        filter(diff == TRUE) %>%
        arrange(qval)
    
    if(nrow(sig_ancom) > 0) {
        cat("  Increased in HFD:", sum(sig_ancom$lfc > 0, na.rm = TRUE), "\n")
        cat("  Decreased in HFD:", sum(sig_ancom$lfc < 0, na.rm = TRUE), "\n\n")
        cat("Top 20 differentially abundant genera:\n")
        print(head(sig_ancom[, c("taxon", "lfc", "qval")], 20))
        
        write.csv(sig_ancom, "ancombc2_diet_effects.csv", row.names = FALSE)
    } else {
        cat("No differentially abundant genera found at alpha = 0.05\n")
        sig_ancom <- data.frame()
    }
} else {
    # If diff column doesn't exist, use q-value threshold
    cat("Using q-value < 0.05 as significance threshold\n")
    
    sig_ancom <- ancom_results %>%
        filter(!is.na(qval) & qval < 0.05) %>%
        arrange(qval)
    
    cat("Significant genera (q < 0.05):", nrow(sig_ancom), "\n")
    if(nrow(sig_ancom) > 0) {
        write.csv(sig_ancom, "ancombc2_diet_effects.csv", row.names = FALSE)
    }
}

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("VISUALIZATION OF ANCOM-BC2 RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# --- Volcano Plot ---
if(nrow(ancom_results) > 0 && all(c("lfc", "pval") %in% colnames(ancom_results))) {
    
    volcano_df <- ancom_results %>%
        mutate(
            neg_log10_p = -log10(pval + 1e-300),
            significant = ifelse(
                "diff" %in% colnames(.), 
                ifelse(diff == TRUE, 
                       ifelse(lfc > 0, "Up in HFD", "Down in HFD"), 
                       "NS"),
                ifelse(!is.na(qval) & qval < 0.05,
                       ifelse(lfc > 0, "Up in HFD", "Down in HFD"),
                       "NS")
            )
        )
    
    p_volcano <- ggplot(volcano_df, aes(x = lfc, y = neg_log10_p)) +
        geom_point(aes(color = significant), alpha = 0.7, size = 2.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
        labs(
            x = "Log Fold Change (HFD vs ND)",
            y = "-log10(P-value)",
            title = "ANCOM-BC2: Differential Abundance (Diet Effect)",
            subtitle = paste0(
                sum(volcano_df$significant == "Up in HFD"), " increased, ",
                sum(volcano_df$significant == "Down in HFD"), " decreased in HFD"
            )
        ) +
        scale_color_manual(
            values = c("Up in HFD" = "#A23B72", 
                       "Down in HFD" = "#2E86AB",
                       "NS" = "grey70"),
            name = ""
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")
    
    # Add labels for top significant taxa
    if(nrow(sig_ancom) > 0) {
        top_label <- sig_ancom %>%
            arrange(qval) %>%
            head(15)
        
        p_volcano <- p_volcano +
            geom_text_repel(
                data = volcano_df %>% filter(taxon %in% top_label$taxon),
                aes(label = taxon),
                size = 3.5,
                max.overlaps = 15,
                box.padding = 0.5,
                fontface = "italic"
            )
    }
    
    ggsave("ancombc2_volcano_diet.pdf", p_volcano, width = 10, height = 8)
    print(p_volcano)
}

# --- Forest Plot ---
if(nrow(sig_ancom) > 0 && all(c("lfc", "se") %in% colnames(sig_ancom))) {
    
    top_forest <- sig_ancom %>%
        arrange(desc(abs(lfc))) %>%
        head(25)
    
    p_forest <- ggplot(top_forest, aes(x = lfc, y = reorder(taxon, lfc))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), 
                       height = 0.2, linewidth = 1, color = "grey40") +
        geom_point(aes(color = lfc > 0), size = 3) +
        labs(
            x = "Log Fold Change (HFD vs ND)",
            y = "",
            title = "ANCOM-BC2: Top 25 Differential Genera",
            subtitle = "Error bars indicate standard error"
        ) +
        scale_color_manual(
            values = c("TRUE" = "#A23B72", "FALSE" = "#2E86AB"),
            labels = c("Decreased in HFD", "Increased in HFD"),
            name = ""
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.y = element_text(face = "italic"),
            legend.position = "bottom"
        )
    
    ggsave("ancombc2_forest_plot.pdf", p_forest, width = 10, height = 8)
    print(p_forest)
}

# --- Heatmap of significant genera ---
if(nrow(sig_ancom) > 0 && exists("genus_rel")) {
    
    sig_genera <- sig_ancom$taxon[sig_ancom$taxon %in% rownames(genus_rel)]
    
    if(length(sig_genera) > 0) {
        sig_genus_rel <- genus_rel[sig_genera, , drop = FALSE]
        sig_genus_log <- log10(sig_genus_rel + 1e-5)
        
        anno_col <- data.frame(
            Diet = metadata_clean$group1,
            Phase = metadata_clean$Phase,
            row.names = colnames(sig_genus_log)
        )
        
        anno_colors <- list(
            Diet = c("ND" = "#2E86AB", "HFD" = "#A23B72"),
            Phase = c("Acute" = "#440154", "Early" = "#3B528B", 
                      "Mid" = "#21908C", "Late" = "#5DC863")
        )
        
        if(length(sig_genera) > 1) {
            pdf("ancombc2_significant_genera_heatmap.pdf", 
                width = 14, height = max(8, length(sig_genera) * 0.35))
            
            pheatmap(sig_genus_log,
                     scale = "row",
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "euclidean",
                     annotation_col = anno_col,
                     annotation_colors = anno_colors,
                     show_rownames = TRUE,
                     fontsize_row = 9,
                     main = paste("ANCOM-BC2 Significant Genera (n =", length(sig_genera), ")"),
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
            dev.off()
            cat("Heatmap saved to ancombc2_significant_genera_heatmap.pdf\n")
        }
    }
}

# ============================================================================
# PHASE-STRATIFIED ANALYSIS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANCOM-BC2: PHASE-STRATIFIED ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

all_phase_sig <- data.frame()
phase_results_list <- list()

for(phase in levels(sample_meta$Phase)) {
    cat("\n", paste(rep("-", 60), collapse = ""), "\n")
    cat("Phase:", phase, "\n")
    
    # Subset TSE
    phase_idx <- colData(tse_genus)$Phase == phase
    tse_phase <- tse_genus[, phase_idx]
    
    # Check group distribution
    groups_present <- table(colData(tse_phase)$group1)
    cat("  Group distribution:", paste(names(groups_present), groups_present, sep = "=", collapse = ", "), "\n")
    
    if(length(unique(colData(tse_phase)$group1)) < 2) {
        cat("  Only one group present. Skipping.\n")
        next
    }
    
    ancom_phase <- tryCatch({
        ancombc2(
            data = tse_phase,
            assay_name = "counts",
            tax_level = "genus",
            fix_formula = "group1 + day_scaled",
            p_adj_method = "holm",
            prv_cut = 0.2,
            lib_cut = 1000,
            group = "group1",
            struc_zero = TRUE,
            neg_lb = TRUE,
            alpha = 0.05,
            n_cl = 4
        )
    }, error = function(e) {
        cat("  Error:", e$message, "\n")
        return(NULL)
    })
    
    if(!is.null(ancom_phase)) {
        phase_results_list[[phase]] <- ancom_phase
        
        res_p <- ancom_phase$res
        lfc_col <- grep("^lfc_group1", colnames(res_p), value = TRUE)
        q_col <- grep("^q_group1", colnames(res_p), value = TRUE)
        diff_col <- grep("^diff_group1", colnames(res_p), value = TRUE)
        
        if(length(diff_col) > 0 && length(lfc_col) > 0) {
            phase_sig <- data.frame(
                taxon = res_p$taxon,
                phase = phase,
                lfc = res_p[[lfc_col]],
                qval = res_p[[q_col]],
                diff = res_p[[diff_col]],
                stringsAsFactors = FALSE
            ) %>%
            filter(diff == TRUE)
        } else {
            # Use q-value threshold
            phase_sig <- data.frame(
                taxon = res_p$taxon,
                phase = phase,
                lfc = res_p[[lfc_col]],
                qval = res_p[[q_col]],
                diff = res_p[[q_col]] < 0.05,
                stringsAsFactors = FALSE
            ) %>%
            filter(diff == TRUE)
        }
        
        cat("  Significant genera:", nrow(phase_sig), "\n")
        
        if(nrow(phase_sig) > 0) {
            all_phase_sig <- rbind(all_phase_sig, phase_sig)
            write.csv(phase_sig, paste0("ancombc2_", phase, ".csv"), row.names = FALSE)
        }
    }
}

# Phase summary
if(nrow(all_phase_sig) > 0) {
    phase_summary <- all_phase_sig %>%
        group_by(phase) %>%
        summarise(
            n_genera = n(),
            up_HFD = sum(lfc > 0),
            down_HFD = sum(lfc < 0),
            .groups = 'drop'
        )
    
    cat("\nPhase-specific summary:\n")
    print(phase_summary)
    write.csv(all_phase_sig, "ancombc2_phase_specific_all.csv", row.names = FALSE)
}

# ============================================================================
# COMPARE WITH MaAsLin3
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("COMPARISON WITH MaAsLin3\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

maaslin_file <- "maaslin3_diet_effect/significant_results.tsv"
if(file.exists(maaslin_file)) {
    maaslin_sig <- read.delim(maaslin_file)
    maaslin_diet <- maaslin_sig %>%
        filter(metadata == "group1") %>%
        mutate(taxon = gsub("^g__", "", feature))
    
    if(nrow(sig_ancom) > 0 && nrow(maaslin_diet) > 0) {
        overlap <- intersect(sig_ancom$taxon, maaslin_diet$taxon)
        ancom_only <- setdiff(sig_ancom$taxon, maaslin_diet$taxon)
        maaslin_only <- setdiff(maaslin_diet$taxon, sig_ancom$taxon)
        
        cat("Method overlap:\n")
        cat("  Both methods:", length(overlap), "\n")
        cat("  ANCOM-BC2 only:", length(ancom_only), "\n")
        cat("  MaAsLin3 only:", length(maaslin_only), "\n")
        
        comparison_df <- data.frame(
            method = c(rep("Both", length(overlap)), 
                       rep("ANCOM-BC2 only", length(ancom_only)),
                       rep("MaAsLin3 only", length(maaslin_only))),
            genus = c(overlap, ancom_only, maaslin_only)
        )
        write.csv(comparison_df, "method_comparison.csv", row.names = FALSE)
        
        if(length(overlap) > 0) {
            cat("\nConcordant genera (detected by both):\n")
            print(overlap)
        }
    }
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SAVING FINAL RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Save all ANCOM-BC2 results
if(nrow(sig_ancom) > 0) {
    write.csv(sig_ancom, "ancombc2_final_results.csv", row.names = FALSE)
}

# Summary report
sink("ancombc2_summary_report.txt")
cat("ANCOM-BC2 ANALYSIS SUMMARY\n")
cat("==========================\n\n")
cat("Date:", date(), "\n\n")
cat("Diet Effect (HFD vs ND):\n")
cat("  Genera tested:", nrow(ancom_results), "\n")
cat("  Significant:", nrow(sig_ancom), "\n")
if(nrow(sig_ancom) > 0) {
    cat("  Up in HFD:", sum(sig_ancom$lfc > 0), "\n")
    cat("  Down in HFD:", sum(sig_ancom$lfc < 0), "\n")
}
cat("\nPhase-stratified results:\n")
if(nrow(all_phase_sig) > 0) {
    print(phase_summary)
}
sink()

# Save workspace
save.image("ancombc2_analysis_complete.RData")

cat("\nAnalysis complete! Generated files:\n")
cat("  - ancombc2_volcano_diet.pdf\n")
if(exists("p_forest")) cat("  - ancombc2_forest_plot.pdf\n")
if(exists("sig_genus_rel")) cat("  - ancombc2_significant_genera_heatmap.pdf\n")
cat("  - ancombc2_diet_effects.csv\n")
cat("  - ancombc2_final_results.csv\n")
if(nrow(all_phase_sig) > 0) cat("  - ancombc2_phase_specific_all.csv\n")
if(exists("comparison_df")) cat("  - method_comparison.csv\n")
cat("  - ancombc2_summary_report.txt\n")
cat("  - ancombc2_analysis_complete.RData\n")