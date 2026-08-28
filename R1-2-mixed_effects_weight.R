#=============================================================================
# R1-2: MIXED-EFFECTS MODEL FOR WEIGHT-ABUNDANCE RELATIONSHIP
#=============================================================================

library(lme4)
library(lmerTest)
library(dplyr)

# Prepare data for the five focal genera
mixed_model_results <- data.frame()

# Define the five genera
five_genera <- c("Lactococcus", "Duncaniella", "Taurinivorans", "Lepagella", "CAG-873")

for(g in five_genera) {
    cat("\nProcessing:", g, "\n")
    
    # Get abundance
    if(g %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[g, ])
    } else if(paste0("g__", g) %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[paste0("g__", g), ])
    } else {
        cat("  Not found. Skipping.\n")
        next
    }
    
    # Create model dataframe
    model_df <- data.frame(
        Abundance = abund,
        LogAbundance = log10(abund + 1e-6),
        Weight = as.numeric(metadata_clean$weight),
        Diet = factor(metadata_clean$group1, levels = c("ND", "HFD")),
        Cage = factor(metadata_clean$cage),
        Day = as.numeric(metadata_clean$Day)
    )
    
    # Add mouse ID if available
    if("mice" %in% colnames(metadata_clean)) {
        model_df$Mouse <- factor(metadata_clean$mice)
    }
    
    # Fit mixed-effects model with fallback options
    model <- NULL
    model_type <- ""
    
    # Attempt 1: Full model with mouse nested in cage
    if("Mouse" %in% colnames(model_df) && is.null(model)) {
        model <- tryCatch({
            m <- lmer(Weight ~ LogAbundance + Diet + Day + (1|Cage/Mouse), 
                      data = model_df, REML = TRUE)
            model_type <- "Cage/Mouse random effects"
            m
        }, error = function(e) NULL)
    }
    
    # Attempt 2: Cage only random effect
    if(is.null(model)) {
        model <- tryCatch({
            m <- lmer(Weight ~ LogAbundance + Diet + Day + (1|Cage), 
                      data = model_df, REML = TRUE)
            model_type <- "Cage random effect"
            m
        }, error = function(e) NULL)
    }
    
    # Attempt 3: Simple linear model (no random effects)
    if(is.null(model)) {
        model <- tryCatch({
            m <- lm(Weight ~ LogAbundance + Diet + Day, data = model_df)
            model_type <- "Simple linear model"
            m
        }, error = function(e) NULL)
    }
    
    if(is.null(model)) {
        cat("  All models failed. Skipping.\n")
        next
    }
    
    cat("  Model used:", model_type, "\n")
    
    # Check for singular fit
    if(inherits(model, "lmerMod")) {
        if(isSingular(model, tol = 1e-4)) {
            cat("  WARNING: Singular fit detected. Simplifying model...\n")
            
            # Remove Day if singular
            model <- tryCatch({
                m <- lmer(Weight ~ LogAbundance + Diet + (1|Cage), 
                          data = model_df, REML = TRUE)
                cat("  Simplified to: Weight ~ LogAbundance + Diet + (1|Cage)\n")
                m
            }, error = function(e) NULL)
        }
    }
    
    if(is.null(model)) {
        cat("  Even simplified model failed. Skipping.\n")
        next
    }
    
    # Extract coefficients safely
    if(inherits(model, "lmerMod")) {
        coef_res <- summary(model)$coefficients
        anova_res <- anova(model)
        
        # Check which terms are present
        terms_present <- rownames(coef_res)
        cat("  Terms in model:", paste(terms_present, collapse = ", "), "\n")
        
        # Extract abundance coefficient
        if("LogAbundance" %in% rownames(coef_res)) {
            abund_coef <- coef_res["LogAbundance", "Estimate"]
            abund_se <- coef_res["LogAbundance", "Std. Error"]
        } else {
            abund_coef <- NA
            abund_se <- NA
        }
        
        # Extract abundance p-value
        if("LogAbundance" %in% rownames(anova_res)) {
            abund_p <- anova_res["LogAbundance", "Pr(>F)"]
        } else {
            abund_p <- NA
        }
        
        # Extract Diet coefficient (HFD vs ND)
        diet_hfd_row <- grep("DietHFD|Diet", rownames(coef_res), value = TRUE)
        if(length(diet_hfd_row) > 0) {
            diet_coef <- coef_res[diet_hfd_row[1], "Estimate"]
            diet_p <- if(diet_hfd_row[1] %in% rownames(anova_res)) 
                anova_res[diet_hfd_row[1], "Pr(>F)"] else NA
        } else {
            diet_coef <- NA
            diet_p <- NA
        }
        
    } else {
        # For lm model
        coef_res <- summary(model)$coefficients
        anova_res <- anova(model)
        
        abund_coef <- if("LogAbundance" %in% rownames(coef_res)) 
            coef_res["LogAbundance", "Estimate"] else NA
        abund_se <- if("LogAbundance" %in% rownames(coef_res))
            coef_res["LogAbundance", "Std. Error"] else NA
        abund_p <- if("LogAbundance" %in% rownames(anova_res))
            anova_res["LogAbundance", "Pr(>F)"] else NA
        
        diet_coef <- if("Diet" %in% rownames(coef_res)) 
            coef_res["Diet", "Estimate"] else NA
        diet_p <- if("Diet" %in% rownames(anova_res))
            anova_res["Diet", "Pr(>F)"] else NA
    }
    
    # Store results
    mixed_model_results <- rbind(mixed_model_results, data.frame(
        Genus = g,
        Model = model_type,
        Abundance_coef = round(abund_coef, 4),
        Abundance_SE = round(abund_se, 4),
        Abundance_p = round(abund_p, 4),
        Diet_coef = round(diet_coef, 4),
        Diet_p = round(diet_p, 4),
        stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  Abundance effect: coef = %.3f, p = %.4f\n",
                abund_coef, abund_p))
    cat(sprintf("  Diet effect: coef = %.3f, p = %.4f\n",
                diet_coef, diet_p))
}

# Add significance labels
mixed_model_results <- mixed_model_results %>%
    mutate(
        Abundance_sig = case_when(
            Abundance_p < 0.001 ~ "***",
            Abundance_p < 0.01 ~ "**",
            Abundance_p < 0.05 ~ "*",
            TRUE ~ "ns"
        ),
        Diet_sig = case_when(
            Diet_p < 0.001 ~ "***",
            Diet_p < 0.01 ~ "**",
            Diet_p < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

cat("\n\n--- Final Mixed-Effects Model Results ---\n")
print(mixed_model_results)
write.csv(mixed_model_results, "Table_mixed_effects_weight.csv", row.names = FALSE)

#=============================================================================
# PER-MOUSE AGGREGATED CORRELATIONS (Simple, Robust Alternative)
#=============================================================================

cat("\n--- Per-mouse aggregated Spearman correlations ---\n")

# Determine mouse identifier
if("mice" %in% colnames(metadata_clean)) {
    mouse_id <- metadata_clean$mice
    cat("Using 'mice' column for mouse ID\n")
} else if("SampleID" %in% colnames(metadata_clean)) {
    # Extract mouse ID from SampleID (e.g., "A51_DAY1" -> "A51")
    mouse_id <- gsub("_DAY.*", "", metadata_clean$SampleID)
    cat("Deriving mouse ID from SampleID\n")
} else {
    mouse_id <- metadata_clean$cage
    cat("Using cage as proxy for mouse\n")
}

per_mouse_df <- data.frame()

for(g in five_genera) {
    # Get abundance
    if(g %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[g, ])
    } else if(paste0("g__", g) %in% rownames(genus_rel)) {
        abund <- as.numeric(genus_rel[paste0("g__", g), ])
    } else {
        next
    }
    
    mouse_data <- data.frame(
        Abundance = abund,
        Weight = as.numeric(metadata_clean$weight),
        Mouse = factor(mouse_id)
    )
    
    # Aggregate per mouse
    mouse_agg <- mouse_data %>%
        group_by(Mouse) %>%
        summarise(
            MeanAbund = mean(Abundance, na.rm = TRUE),
            MeanWeight = mean(Weight, na.rm = TRUE),
            .groups = 'drop'
        )
    
    if(nrow(mouse_agg) >= 5) {
        ct <- cor.test(mouse_agg$MeanAbund, mouse_agg$MeanWeight, 
                       method = "spearman", exact = FALSE)
        
        per_mouse_df <- rbind(per_mouse_df, data.frame(
            Genus = g,
            N_mice = nrow(mouse_agg),
            Spearman_rho = round(ct$estimate, 4),
            P_value = round(ct$p.value, 4)
        ))
    }
}

# Adjust p-values
if(nrow(per_mouse_df) > 0) {
    per_mouse_df$P_adj <- round(p.adjust(per_mouse_df$P_value, method = "BH"), 4)
    per_mouse_df$Significance <- case_when(
        per_mouse_df$P_adj < 0.001 ~ "***",
        per_mouse_df$P_adj < 0.01 ~ "**",
        per_mouse_df$P_adj < 0.05 ~ "*",
        TRUE ~ "ns"
    )
}

cat("\nPer-mouse aggregated correlations:\n")
print(per_mouse_df)
write.csv(per_mouse_df, "Table_per_mouse_correlations.csv", row.names = FALSE)

cat("\nDone!\n")