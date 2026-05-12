#=============================================================================
# IMPROVED CO-OCCURRENCE NETWORK VISUALIZATION
# Clear labeling of five focal genera, no overlapping nodes
#=============================================================================


library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FIXED NETWORK - CHECKING NAME MATCHING\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# CHECK NAMES
# ============================================================================

# First, check what names are actually in the network
cat("Five focal genera:\n")
cat("  ", paste(five_genera, collapse = ", "), "\n\n")

cat("All node names in network:\n")
cat("  ", paste(V(g)$name, collapse = ", "), "\n\n")

cat("Cleaned node names:\n")
cat("  ", paste(V(g)$name_clean, collapse = ", "), "\n\n")

# Check which focal genera match
for(fg in five_genera) {
    # Try exact match with cleaned names
    match_clean <- which(V(g)$name_clean == fg)
    
    # Try matching original names
    match_orig <- which(V(g)$name == fg)
    
    # Try with g__ prefix
    match_g <- which(V(g)$name == paste0("g__", fg))
    
    # Try grep
    match_grep <- grep(fg, V(g)$name, ignore.case = TRUE)
    
    cat(sprintf("Matching '%s':\n", fg))
    cat(sprintf("  Exact (clean): %s\n", if(length(match_clean) > 0) V(g)$name[match_clean] else "NOT FOUND"))
    cat(sprintf("  Exact (orig):  %s\n", if(length(match_orig) > 0) V(g)$name[match_orig] else "NOT FOUND"))
    cat(sprintf("  g__ prefix:    %s\n", if(length(match_g) > 0) V(g)$name[match_g] else "NOT FOUND"))
    cat(sprintf("  Grep:          %s\n", if(length(match_grep) > 0) paste(V(g)$name[match_grep], collapse=", ") else "NOT FOUND"))
    cat("\n")
}

# ============================================================================
# FIX MATCHING LOGIC
# ============================================================================

# Try multiple matching strategies
match_focal <- function(node_names, focal_list) {
    matched <- rep(FALSE, length(node_names))
    
    for(fg in focal_list) {
        # Strategy 1: exact match (cleaned)
        idx <- which(node_names == fg)
        if(length(idx) > 0) { matched[idx] <- TRUE; next }
        
        # Strategy 2: match with g__ prefix
        idx <- which(node_names == paste0("g__", fg))
        if(length(idx) > 0) { matched[idx] <- TRUE; next }
        
        # Strategy 3: case-insensitive grep
        idx <- grep(fg, node_names, ignore.case = TRUE)
        if(length(idx) > 0) { matched[idx[1]] <- TRUE; next }
    }
    
    return(matched)
}

# Apply matching
V(g)$is_focal <- match_focal(V(g)$name_clean, five_genera) | 
    match_focal(V(g)$name, five_genera)

cat("Focal genera found after fixing:", sum(V(g)$is_focal), "\n")
cat("  ", paste(V(g)$name_clean[V(g)$is_focal], collapse = ", "), "\n")

# If still 0, manually match
if(sum(V(g)$is_focal) == 0) {
    cat("\nManual matching required. Checking partial matches...\n")
    
    # Check each focal genus against all node names
    for(fg in five_genera) {
        # Remove special characters and try matching
        fg_simple <- gsub("[^a-zA-Z0-9]", "", fg)
        
        for(i in 1:length(V(g)$name)) {
            node_simple <- gsub("[^a-zA-Z0-9]", "", V(g)$name[i])
            
            if(tolower(fg_simple) == tolower(node_simple) ||
               grepl(tolower(fg_simple), tolower(node_simple)) ||
               grepl(tolower(node_simple), tolower(fg_simple))) {
                cat(sprintf("  Match: '%s' <-> '%s'\n", fg, V(g)$name[i]))
                V(g)$is_focal[i] <- TRUE
            }
        }
    }
    
    cat("\nFocal genera found after manual matching:", sum(V(g)$is_focal), "\n")
}

# ============================================================================
# PROCEED WITH IDENTIFIED FOCAL GENERA
# ============================================================================

if(sum(V(g)$is_focal) == 0) {
    cat("\nERROR: Still cannot find focal genera. Using top 5 by degree instead.\n")
    # Fallback: use top 5 by degree
    top5_idx <- order(V(g)$degree, decreasing = TRUE)[1:5]
    V(g)$is_focal <- FALSE
    V(g)$is_focal[top5_idx] <- TRUE
    cat("Using as focal:", paste(V(g)$name_clean[top5_idx], collapse = ", "), "\n")
}

# Update properties
V(g)$label <- ifelse(V(g)$is_focal, V(g)$name_clean, "")
V(g)$category <- ifelse(V(g)$is_focal, "Focal",
                        ifelse(V(g)$degree >= quantile(V(g)$degree, 0.85), 
                               "Hub", "Other"))

focal_node_names <- V(g)$name[V(g)$is_focal]

# ============================================================================
# CUSTOM LAYOUT
# ============================================================================

set.seed(123)
layout_use <- layout_with_fr(g, niter = 2000)

focal_idx <- which(V(g)$is_focal)
if(length(focal_idx) > 0) {
    n_focal <- length(focal_idx)
    angles <- seq(0, 2 * pi, length.out = n_focal + 1)[1:n_focal]
    radius <- max(dist(layout_use)) * 0.3
    
    for(i in 1:n_focal) {
        layout_use[focal_idx[i], 1] <- radius * cos(angles[i])
        layout_use[focal_idx[i], 2] <- radius * sin(angles[i])
    }
}

# ============================================================================
# CREATE DATA FRAMES
# ============================================================================

nodes_df <- data.frame(
    x = layout_use[, 1],
    y = layout_use[, 2],
    name = V(g)$name_clean,
    degree = V(g)$degree,
    is_focal = V(g)$is_focal,
    category = V(g)$category,
    label = V(g)$label,
    stringsAsFactors = FALSE
)

edges_df <- as_data_frame(g, what = "edges")
edges_df$x_start <- layout_use[match(edges_df$from, V(g)$name), 1]
edges_df$y_start <- layout_use[match(edges_df$from, V(g)$name), 2]
edges_df$x_end <- layout_use[match(edges_df$to, V(g)$name), 1]
edges_df$y_end <- layout_use[match(edges_df$to, V(g)$name), 2]
edges_df$direction <- ifelse(edges_df$Correlation > 0, "Positive", "Negative")

# ============================================================================
# PLOT
# ============================================================================

p1 <- ggplot() +
    # Edges
    geom_segment(
        data = edges_df,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = direction, alpha = abs(Correlation)),
        linewidth = 0.3
    ) +
    scale_color_manual(
        values = c("Positive" = "#E8A0A0", "Negative" = "#A0C4E8"),
        name = "Correlation"
    ) +
    scale_alpha_continuous(range = c(0.1, 0.6), guide = "none") +
    
    # Non-focal nodes
    geom_point(
        data = nodes_df %>% filter(!is_focal),
        aes(x = x, y = y, size = degree, fill = category),
        shape = 21, color = "white", stroke = 0.3, alpha = 0.7
    ) +
    # Focal nodes (on top, larger border)
    geom_point(
        data = nodes_df %>% filter(is_focal),
        aes(x = x, y = y, size = degree),
        shape = 21, fill = "#E63946", color = "white", stroke = 1.5, alpha = 1
    ) +
    scale_fill_manual(
        values = c("Focal" = "#E63946", "Hub" = "#F4A261", "Other" = "#ADB5BD"),
        guide = "none"
    ) +
    scale_size_continuous(
        range = c(2, 15),
        name = "Degree (connections)",
        guide = guide_legend(override.aes = list(fill = "grey60", shape = 21))
    ) +
    
    # White background for focal labels
    ggrepel::geom_label_repel(
        data = nodes_df %>% filter(is_focal),
        aes(x = x, y = y, label = label),
        fill = "white",
        alpha = 0.8,
        color = "#8B0000",
        fontface = "bold.italic",
        size = 5,
        label.size = NA,
        label.padding = unit(0.2, "lines"),
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = NA,
        show.legend = FALSE
    ) +
    
    labs(
        title = "Co-occurrence Network of Gut Microbiota",
        subtitle = paste0(
            "Spearman |r| > 0.5, p.adj < 0.05 | ",
            nrow(nodes_df), " genera, ", nrow(edges_df), " edges"
        ),
        caption = "Red nodes: focal genera | Orange: hubs | Grey: other | Red edges: positive | Blue: negative"
    ) +
    theme_void() +
    theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5),
        plot.caption = element_text(size = 9, color = "grey60", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9)
    ) +
    coord_fixed()

ggsave("network_focal_labeled.pdf", p1, width = 14, height = 12, dpi = 300)
print(p1)

# ============================================================================
# CONNECTIONS SUMMARY
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FOCAL GENERA CONNECTIONS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

for(fn in focal_node_names) {
    fn_clean <- clean_genus_name(fn)
    neigh <- names(neighbors(g, fn))
    neigh_clean <- clean_genus_name(neigh)
    
    cat(sprintf("\n%s (Degree: %d):\n", fn_clean, length(neigh)))
    cat("  ", paste(neigh_clean, collapse = ", "), "\n")
}

save(p1, g, nodes_df, edges_df, file = "network_final.RData")
cat("\nDone! Check network_focal_labeled.pdf\n")