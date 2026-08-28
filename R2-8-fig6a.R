#=============================================================================
# R2-8: REVISED COMMUNITY STATE PLOT WITH TRANSITIONS
#=============================================================================


library(ggplot2)
library(dplyr)
library(cluster)
library(patchwork)

cat("Running PAM clustering and creating community state plot...\n\n")

#=============================================================================
# 1. RUN PAM CLUSTERING
#=============================================================================

# Check if bray_dist exists
if(!exists("bray_dist")) {
    cat("Calculating Bray-Curtis distance...\n")
    bray_dist <- vegdist(t(genus_rel), method = "bray")
}

cat("Bray-Curtis distance matrix:", nrow(bray_dist), "x", ncol(bray_dist), "\n")

# Run PAM clustering with k=2
set.seed(123)
pam_final <- pam(bray_dist, k = 2, diss = TRUE)

cat("PAM clustering complete\n")
cat("Number of clusters:", length(unique(pam_final$clustering)), "\n")
cat("Silhouette width:", round(pam_final$silinfo$avg.width, 4), "\n")
cat("Cluster assignments:\n")
print(table(pam_final$clustering))

#=============================================================================
# 2. BUILD STATE DATA
#=============================================================================

state_data <- data.frame(
    SampleID = colnames(genus_rel),
    Day = as.numeric(metadata_clean$Day),
    Diet = factor(metadata_clean$group1, levels = c("ND", "HFD")),
    Cage = factor(metadata_clean$cage),
    CommunityState = as.factor(pam_final$clustering),
    stringsAsFactors = FALSE
)

cat("\nState data built:\n")
cat("  Rows:", nrow(state_data), "\n")
cat("  Columns:", paste(colnames(state_data), collapse = ", "), "\n")
cat("  States:", paste(levels(state_data$CommunityState), collapse = ", "), "\n")
cat("  Cages:", paste(levels(state_data$Cage), collapse = ", "), "\n")

# Show state distribution by diet
cat("\nState distribution by diet:\n")
print(table(state_data$Diet, state_data$CommunityState))

#=============================================================================
# 3. CALCULATE TRANSITIONS
#=============================================================================

transition_df <- data.frame()

for(cage_id in levels(state_data$Cage)) {
    cage_data <- state_data %>%
        filter(Cage == cage_id) %>%
        arrange(Day)
    
    if(nrow(cage_data) > 1) {
        for(i in 1:(nrow(cage_data)-1)) {
            transition_df <- rbind(transition_df, data.frame(
                Cage = cage_id,
                Diet = cage_data$Diet[i],
                Day = cage_data$Day[i],
                NextDay = cage_data$Day[i+1],
                From = as.character(cage_data$CommunityState[i]),
                To = as.character(cage_data$CommunityState[i+1]),
                Changed = cage_data$CommunityState[i] != cage_data$CommunityState[i+1],
                stringsAsFactors = FALSE
            ))
        }
    }
}

cat("\nTotal transition pairs:", nrow(transition_df), "\n")
cat("State changes:", sum(transition_df$Changed), "\n")

#=============================================================================
# 4. TRANSITION RATES
#=============================================================================

transition_rates <- transition_df %>%
    group_by(Diet) %>%
    summarise(
        N_transitions = sum(Changed),
        N_total = n(),
        Rate = N_transitions / N_total * 100,
        .groups = 'drop'
    )

cat("\nTransition rates by diet:\n")
print(transition_rates)

#=============================================================================
# 5. STATE DYNAMICS PLOT
#=============================================================================

state_colors <- c("1" = "#FF6B6B", "2" = "#4ECDC4")

p_states <- ggplot() +
    # State tiles
    geom_tile(
        data = state_data,
        aes(x = Day, y = Cage, fill = CommunityState),
        color = "white", linewidth = 1
    ) +
    # Transition arrows
    geom_segment(
        data = transition_df %>% filter(Changed),
        aes(x = Day + 0.5, xend = NextDay - 0.5,
            y = Cage, yend = Cage),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        color = "black", linewidth = 0.8, alpha = 0.8
    ) +
    # Facet by Diet
    facet_wrap(~ Diet, scales = "free_y", ncol = 1) +
    # Scales
    scale_fill_manual(
        values = state_colors,
        name = "Community\nState",
        labels = c("State 1", "State 2")
    ) +
    scale_x_continuous(
        breaks = sort(unique(state_data$Day)),
        labels = paste0("D", sort(unique(state_data$Day)))
    ) +
    # Labels
    labs(
        title = "Community State Dynamics During HFD-Induced Obesity",
        subtitle = paste0(
            "PAM clustering (k=2) of Bray-Curtis distances | Silhouette = ",
            round(pam_final$silinfo$avg.width, 3), "\n",
            "Arrows = state transitions between consecutive timepoints"
        ),
        x = "Day",
        y = "Cage"
    ) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, color = "grey40", lineheight = 1.2),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid = element_blank(),
        legend.position = "right"
    )

ggsave("Revised_Figure6A_state_dynamics.pdf", p_states, 
       width = 12, height = 7, dpi = 300)
print(p_states)

#=============================================================================
# 6. TRANSITION FREQUENCY PLOT
#=============================================================================

p_freq <- ggplot(transition_rates, aes(x = Diet, y = Rate, fill = Diet)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.5) +
    geom_text(aes(label = paste0(round(Rate, 1), "%")), 
              vjust = -0.5, size = 6, fontface = "bold") +
    scale_fill_manual(values = c("ND" = "#2E86AB", "HFD" = "#A23B72"), 
                      guide = "none") +
    labs(
        title = "State Transition Frequency",
        subtitle = paste0(
            "HFD: ", transition_rates$N_transitions[transition_rates$Diet == "HFD"], 
            " transitions / ", transition_rates$N_total[transition_rates$Diet == "HFD"], " pairs\n",
            "ND: ", transition_rates$N_transitions[transition_rates$Diet == "ND"], 
            " transitions / ", transition_rates$N_total[transition_rates$Diet == "ND"], " pairs"
        ),
        x = "",
        y = "Transition Rate (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40", lineheight = 1.3)
    ) +
    ylim(0, max(transition_rates$Rate) * 1.25)

ggsave("Revised_Figure6A_transition_frequency.pdf", p_freq, 
       width = 5, height = 6, dpi = 300)
print(p_freq)

#=============================================================================
# 7. COMBINED FIGURE
#=============================================================================

p_combined <- p_states / p_freq +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
        title = "Community State Dynamics and Transition Frequency",
        tag_levels = "A",
        theme = theme(plot.title = element_text(face = "bold", size = 16))
    )

ggsave("Figure6A_combined.pdf", p_combined, width = 12, height = 10, dpi = 300)
print(p_combined)

#=============================================================================
# 8. SAVE
#=============================================================================

save(pam_final, state_data, transition_df, transition_rates,
     p_states, p_freq, p_combined,
     file = "revised_figure6A_data.RData")

write.csv(transition_df, "Table_state_transitions.csv", row.names = FALSE)
write.csv(transition_rates, "Table_transition_rates.csv", row.names = FALSE)

cat("\nDone! Files saved:\n")
cat("  - Revised_Figure6A_state_dynamics.pdf\n")
cat("  - Revised_Figure6A_transition_frequency.pdf\n")
cat("  - Figure6A_combined.pdf\n")
cat("  - Table_state_transitions.csv\n")
cat("  - Table_transition_rates.csv\n")