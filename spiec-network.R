#=============================================================================
# R1-3a & R2-6: SPIEC-EASI NETWORK VISUALIZATION
#=============================================================================

library(SpiecEasi)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(ggrepel)

cat("Running SPIEC-EASI network analysis...\n")

#=============================================================================
# 1. RUN SPIEC-EASI
#=============================================================================

# Prepare data (samples × taxa)
se_data <- t(genus_counts)
se_data <- round(se_data)
mode(se_data) <- "integer"

# Check data
cat("Data dimensions:", dim(se_data), "\n")

# Run SPIEC-EASI with MB method
set.seed(123)
se_mb <- tryCatch({
    spiec.easi(
        se_data,
        method = "mb",
        lambda.min.ratio = 1e-3,
        nlambda = 50,
        pulsar.params = list(rep.num = 30, thresh = 0.05)
    )
}, error = function(e) {
    cat("MB method failed:", e$message, "\n")
    cat("Trying glasso method...\n")
    spiec.easi(
        se_data,
        method = "glasso",
        lambda.min.ratio = 1e-3,
        nlambda = 50,
        pulsar.params = list(rep.num = 30, thresh = 0.05)
    )
})

cat("SPIEC-EASI completed. Class:", class(se_mb), "\n")

#=============================================================================
# 2. EXTRACT NETWORK (FIXED)
#=============================================================================

# The output of spiec.easi is a pulsar object
# Need to use getRefit() for the final network

# Get refit network
se_refit <- getRefit(se_mb)

cat("Refit object class:", class(se_refit), "\n")
cat("Refit object names:", paste(names(se_refit), collapse = ", "), "\n")

# Extract adjacency matrix
# For MB method, the network is in se_refit$est$path or use getOptMerge
se_adj <- tryCatch({
    # Try getOptMerge first
    net <- getOptMerge(se_mb)
    as.matrix(net)
}, error = function(e1) {
    # Try from refit
    tryCatch({
        as.matrix(se_refit$est$path[[se_refit$est$select]])
    }, error = function(e2) {
        # Try direct refit
        tryCatch({
            as.matrix(se_refit)
        }, error = function(e3) {
            cat("All extraction methods failed. Trying alternate approach...\n")
            # Manual extraction from pulsar
            opt_idx <- se_mb$est$select
            as.matrix(se_mb$est$path[[opt_idx]])
        })
    })
})

# Set row/column names
genus_names <- colnames(se_data)
rownames(se_adj) <- genus_names
colnames(se_adj) <- genus_names

cat("Adjacency matrix dimensions:", dim(se_adj), "\n")
cat("Non-zero edges:", sum(se_adj != 0), "\n")

#=============================================================================
# 3. CALCULATE DEGREE CENTRALITY
#=============================================================================

se_degree <- rowSums(abs(se_adj) > 0)

se_hub_df <- data.frame(
    Genus = names(se_degree),
    Degree = as.numeric(se_degree),
    stringsAsFactors = FALSE
) %>%
    arrange(desc(Degree)) %>%
    mutate(
        Genus_clean = gsub("^g__", "", Genus),
        Genus_clean = gsub("_", " ", Genus_clean),
        Rank = 1:n(),
        is_focal = Genus_clean %in% five_genera | 
            gsub("^g__", "", Genus) %in% five_genera
    )

cat("\nSPIEC-EASI hub ranking (top 15):\n")
print(head(se_hub_df, 15))

cat("\nFocal genera in SPIEC-EASI network:\n")
print(se_hub_df %>% filter(is_focal))

#=============================================================================
# 4. CREATE EDGE LIST FOR PLOTTING
#=============================================================================

edge_list <- data.frame()
for(i in 1:(nrow(se_adj)-1)) {
    for(j in (i+1):ncol(se_adj)) {
        if(se_adj[i, j] != 0) {
            edge_list <- rbind(edge_list, data.frame(
                From = rownames(se_adj)[i],
                To = colnames(se_adj)[j],
                Weight = se_adj[i, j]
            ))
        }
    }
}

cat("\nNetwork edges:", nrow(edge_list), "\n")
cat("Positive edges:", sum(edge_list$Weight > 0), "\n")
cat("Negative edges:", sum(edge_list$Weight < 0), "\n")

#=============================================================================
# 5. CREATE NETWORK PLOT (Using ggplot2 directly if ggraph fails)
#=============================================================================

# Create igraph object
if(nrow(edge_list) > 0) {
    g_se <- graph_from_data_frame(edge_list, directed = FALSE,
                                   vertices = se_hub_df$Genus)
    
    V(g_se)$name_clean <- gsub("^g__", "", V(g_se)$name)
    V(g_se)$name_clean <- gsub("_", " ", V(g_se)$name_clean)
    V(g_se)$degree <- degree(g_se)
    V(g_se)$is_focal <- V(g_se)$name_clean %in% five_genera
    V(g_se)$category <- ifelse(V(g_se)$is_focal, "Focal",
                               ifelse(V(g_se)$degree >= quantile(V(g_se)$degree, 0.85, na.rm = TRUE),
                                      "Hub", "Other"))
    
    # Layout
    set.seed(123)
    layout_se <- layout_with_fr(g_se, niter = 2000)
    
    # Spread focal nodes
    focal_idx <- which(V(g_se)$is_focal)
    if(length(focal_idx) > 0) {
        n_focal <- length(focal_idx)
        angles <- seq(0, 2 * pi, length.out = n_focal + 1)[1:n_focal]
        radius <- max(dist(layout_se)) * 0.3
        
        for(i in 1:n_focal) {
            layout_se[focal_idx[i], 1] <- radius * cos(angles[i])
            layout_se[focal_idx[i], 2] <- radius * sin(angles[i])
        }
    }
    
    # Prepare node dataframe
    nodes_df_se <- data.frame(
        x = layout_se[, 1],
        y = layout_se[, 2],
        name = V(g_se)$name_clean,
        degree = V(g_se)$degree,
        is_focal = V(g_se)$is_focal,
        category = V(g_se)$category,
        stringsAsFactors = FALSE
    )
    
    # Prepare edge dataframe
    edges_df_se <- data.frame(
        x_start = layout_se[match(edge_list$From, V(g_se)$name), 1],
        y_start = layout_se[match(edge_list$From, V(g_se)$name), 2],
        x_end = layout_se[match(edge_list$To, V(g_se)$name), 1],
        y_end = layout_se[match(edge_list$To, V(g_se)$name), 2],
        Weight = edge_list$Weight,
        Direction = ifelse(edge_list$Weight > 0, "Positive", "Negative"),
        stringsAsFactors = FALSE
    )
    
    #=========================================================================
    # PLOT: Main SPIEC-EASI Network
    #=========================================================================
    
    p_spieceasi <- ggplot() +
        # Edges
        geom_segment(data = edges_df_se,
                     aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                         color = Direction, alpha = abs(Weight)),
                     linewidth = 0.4) +
        scale_color_manual(values = c("Positive" = "#E74C3C", "Negative" = "#3498DB"),
                           name = "Association") +
        scale_alpha_continuous(range = c(0.1, 0.7), guide = "none") +
        
        # Non-focal nodes
        geom_point(data = nodes_df_se %>% filter(!is_focal),
                   aes(x = x, y = y, size = degree, fill = category),
                   shape = 21, color = "white", stroke = 0.3, alpha = 0.7) +
        # Focal nodes
        geom_point(data = nodes_df_se %>% filter(is_focal),
                   aes(x = x, y = y, size = degree),
                   shape = 21, fill = "#E63946", color = "white", 
                   stroke = 1.5, alpha = 1) +
        scale_fill_manual(values = c("Focal" = "#E63946", "Hub" = "#F4A261", "Other" = "#90A4AE"),
                          guide = "none") +
        scale_size_continuous(range = c(2, 14), name = "Degree") +
        
        # Labels for focal genera
        geom_label_repel(data = nodes_df_se %>% filter(is_focal),
                         aes(x = x, y = y, label = name),
                         fill = "white", alpha = 0.9,
                         color = "#8B0000", fontface = "bold.italic",
                         size = 4.5, box.padding = 0.8,
                         point.padding = 0.5, segment.color = "#E63946",
                         max.overlaps = 10, show.legend = FALSE) +
        
        labs(
            title = "SPIEC-EASI Compositionally Aware Network",
            subtitle = paste0(
                "MB method | ", nrow(nodes_df_se), " nodes | ", nrow(edges_df_se), " edges\n",
                "Red nodes: focal genera"
            ),
            caption = "Red edges: positive partial correlations | Blue edges: negative partial correlations"
        ) +
        theme_void() +
        theme(
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
            plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5, lineheight = 1.2),
            plot.caption = element_text(size = 9, color = "grey50", hjust = 0.5),
            legend.position = "right"
        ) +
        coord_fixed()
    
    ggsave("Supplementary_Figure_spieceasi_network.pdf", p_spieceasi,
           width = 12, height = 10, dpi = 300)
    print(p_spieceasi)
    
} else {
    cat("\nNo edges found. The network may be too sparse.\n")
    cat("Try adjusting lambda.min.ratio to a smaller value.\n")
}

#=============================================================================
# 6. COMPARISON: SPEARMAN vs SPIEC-EASI
#=============================================================================

# Calculate Spearman degrees for comparison
spearman_cor <- cor(t(genus_rel), method = "spearman", use = "pairwise.complete.obs")
spearman_adj <- abs(spearman_cor) > 0.5
diag(spearman_adj) <- 0
spearman_degree <- rowSums(spearman_adj)

spearman_hubs <- data.frame(
    Genus = names(spearman_degree),
    Spearman_Degree = as.numeric(spearman_degree),
    stringsAsFactors = FALSE
) %>%
    mutate(Genus_clean = gsub("^g__", "", Genus))

# Merge
comparison_df <- se_hub_df %>%
    select(Genus_clean, SE_Degree = Degree, SE_Rank = Rank, is_focal) %>%
    left_join(spearman_hubs %>% select(Genus_clean, Spearman_Degree),
              by = "Genus_clean") %>%
    mutate(Spearman_Rank = rank(-Spearman_Degree, ties.method = "first"))

cat("\nComparison: Spearman vs SPIEC-EASI (focal genera):\n")
print(comparison_df %>% filter(is_focal))

# Scatter plot
p_compare <- ggplot(comparison_df, aes(x = Spearman_Degree, y = SE_Degree)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = is_focal, size = is_focal), alpha = 0.7) +
    geom_text_repel(data = comparison_df %>% filter(is_focal),
                    aes(label = Genus_clean),
                    fontface = "italic", size = 4, color = "#E63946",
                    box.padding = 0.7) +
    scale_color_manual(values = c("TRUE" = "#E63946", "FALSE" = "grey60"),
                       guide = "none") +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2), guide = "none") +
    labs(title = "Hub Degree: Spearman vs SPIEC-EASI",
         subtitle = "Focal genera shown in red",
         x = "Spearman Degree",
         y = "SPIEC-EASI Degree") +
    theme_bw(base_size = 12)

ggsave("Supplementary_Figure_network_comparison.pdf", p_compare,
       width = 8, height = 7)
print(p_compare)

#=============================================================================
# 7. SAVE
#=============================================================================

write.csv(se_hub_df, "Table_spieceasi_hubs.csv", row.names = FALSE)
write.csv(comparison_df, "Table_network_comparison.csv", row.names = FALSE)

save(se_mb, se_adj, se_hub_df, comparison_df, 
     file = "spiec_easi_network_results.RData")

cat("\nSPIEC-EASI analysis complete!\n")