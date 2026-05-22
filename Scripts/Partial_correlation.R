library(tidyverse)
library(openxlsx)
library(ppcor)

# Set working and output directories
output_path <- "D:/output_folder"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
setwd(output_path)

# Create output folders
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Load data
data <- read.csv(
  "D:/input_folder/data.csv",
  stringsAsFactors = FALSE
)


cat("Number of rows in the original data:", nrow(data), "\n")
cat("Number of columns in the original data:", ncol(data), "\n")

# Define brain-region list
brain_regions <- c(
  "Quant_lh_Hippocampus", "Quant_rh_Hippocampus",
  "Quant_lh_Amygdala", "Quant_rh_Amygdala",
  "Quant_lh_Accumbens", "Quant_rh_Accumbens",
  "Quant_lh_Anteriorcingulate", "Quant_rh_Anteriorcingulate",
  "Quant_lh_Posteriorcingulate", "Quant_rh_Posteriorcingulate",
  "Quant_lh_Entorhinal", "Quant_rh_Entorhinal",
  "Quant_lh_Orbitofrontal", "Quant_rh_Orbitofrontal"
)

# Define simplified brain-region names
brain_regions_simple <- c(
  "Hippocampus_L", "Hippocampus_R", "Amygdala_L", "Amygdala_R", 
  "Accumbens_L", "Accumbens_R", "Anteriorcingulate_L", "Anteriorcingulate_R",
  "Posteriorcingulate_L", "Posteriorcingulate_R", "Entorhinal_L", "Entorhinal_R",
  "Orbitofrontal_L", "Orbitofrontal_R"
)

# Define clinical scales
all_scales <- c("MOCA", "SDMT", "CVLT", "BVMT", "PASAT", "EDSS", "CASE")

# Define disease groups available for each scale
scale_diseases <- list(
  "MOCA" = c("MS", "AQP4Pos_NMOSD", "MOGAD","NMDA", "LGI1", "GAD65"),
  "SDMT" = c("MS", "AQP4Pos_NMOSD"),
  "CVLT" = c("MS", "AQP4Pos_NMOSD"),
  "BVMT" = c("MS", "AQP4Pos_NMOSD"),
  "PASAT" = c("MS", "AQP4Pos_NMOSD"),
  "EDSS" = c("MS", "AQP4Pos_NMOSD", "MOGAD"),
  "CASE" = c("NMDA", "LGI1", "GAD65")
)


all_diseases <- unique(unlist(scale_diseases))


data_filtered <- data %>%
  filter(Diagnosis %in% all_diseases) %>%
  mutate(
    across(all_of(all_scales), ~ ifelse(. %in% c("", "NA", "N/A", NA), NA, as.numeric(as.character(.)))),
    Age = as.numeric(as.character(Age)),
    Sex = as.character(Sex),
    Sex_numeric = case_when(
      tolower(Sex) %in% c("m", "male", "1", "male ", " male") ~ 1,
      tolower(Sex) %in% c("f", "female", "0", "female ", " female") ~ 0,
      TRUE ~ NA_real_
    ),
    across(all_of(brain_regions), as.numeric)
  )


sample_counts_list <- list()
for (scale in all_scales) {
  scale_counts <- data_filtered %>%
    filter(Diagnosis %in% scale_diseases[[scale]]) %>%
    group_by(Diagnosis) %>%
    summarise(
      n = n(),
      n_valid = sum(!is.na(.data[[scale]]) & !is.na(Age) & !is.na(Sex_numeric))
    ) %>%
    mutate(Scale = scale)
  sample_counts_list[[scale]] <- scale_counts
}
sample_counts_all <- bind_rows(sample_counts_list)
print("Sample counts by disease group after requiring age and sex information:")
print(sample_counts_all)

# 2. Partial correlation function with FDR correction across brain regions
perform_partial_correlation_by_scale <- function() {
  all_results <- data.frame()
  
  for (scale_name in all_scales) {
    diseases_for_scale <- scale_diseases[[scale_name]]
    
    for (disease in diseases_for_scale) {
      disease_data <- data_filtered %>%
        filter(Diagnosis == disease)
      
      scale_results <- data.frame()
      
      for (i in seq_along(brain_regions)) {
        region <- brain_regions[i]
        region_simple <- brain_regions_simple[i]
        
        cols_needed <- c(scale_name, region, "Age", "Sex_numeric")
        missing_cols <- setdiff(cols_needed, names(disease_data))
        if (length(missing_cols) > 0) {
          cat(paste("Warning: missing columns in", disease, ":", paste(missing_cols, collapse = ", "), "\n"))
          next
        }
        
        temp_data <- disease_data[, cols_needed]
        complete_data <- na.omit(temp_data)
        n_complete <- nrow(complete_data)
        
        if (n_complete >= 5) {
          scale_values <- complete_data[[scale_name]]
          region_values <- complete_data[[region]]
          age_values <- complete_data[["Age"]]
          sex_values <- complete_data[["Sex_numeric"]]
          
          tryCatch({
            pcor_result <- pcor.test(
              x = scale_values,
              y = region_values,
              z = data.frame(age_values, sex_values),
              method = "spearman"
            )
            
            temp_result <- data.frame(
              Disease = disease,
              Scale = scale_name,
              BrainRegion = region,
              BrainRegionSimple = region_simple,
              n = n_complete,
              r_value = pcor_result$estimate,
              p_value = pcor_result$p.value,
              stringsAsFactors = FALSE
            )
          }, error = function(e) {
            cat(paste("Error in", disease, "-", scale_name, "-", region, ":", e$message, "\n"))
            temp_result <- data.frame(
              Disease = disease,
              Scale = scale_name,
              BrainRegion = region,
              BrainRegionSimple = region_simple,
              n = n_complete,
              r_value = NA,
              p_value = NA,
              stringsAsFactors = FALSE
            )
          })
        } else {
          temp_result <- data.frame(
            Disease = disease,
            Scale = scale_name,
            BrainRegion = region,
            BrainRegionSimple = region_simple,
            n = n_complete,
            r_value = NA,
            p_value = NA,
            stringsAsFactors = FALSE
          )
        }
        scale_results <- rbind(scale_results, temp_result)
      }
      
      scale_results$p_adj <- p.adjust(scale_results$p_value, method = "fdr")
      all_results <- rbind(all_results, scale_results)
    }
  }
  return(all_results)
}

# 3. Partial correlation analysis
cat("\nStarting partial correlation analysis controlling for age and sex, with FDR correction...\n")
all_results <- perform_partial_correlation_by_scale()

results_by_scale <- list()
for (scale in all_scales) {
  results_by_scale[[scale]] <- all_results %>% filter(Scale == scale)
  cat(paste0(scale, " partial correlation analysis completed: "), nrow(results_by_scale[[scale]]), " correlations tested\n")
}

# 4. Add significance markers
add_significance_stars <- function(p_adj) {
  stars <- rep("", length(p_adj))
  stars[!is.na(p_adj) & p_adj < 0.05] <- "*"
  stars[!is.na(p_adj) & p_adj < 0.01] <- "**"
  stars[!is.na(p_adj) & p_adj < 0.001] <- "***"
  return(stars)
}

# 5. Create display text
create_display_text <- function(results) {
  results$r_text <- ifelse(is.na(results$r_value), "NA", sprintf("%.2f", results$r_value))
  results$sig_text <- results$significance
  return(results)
}

for (scale in all_scales) {
  if (nrow(results_by_scale[[scale]]) > 0) {
    results_by_scale[[scale]]$significance <- add_significance_stars(results_by_scale[[scale]]$p_adj)
    results_by_scale[[scale]] <- create_display_text(results_by_scale[[scale]])
  }
}

# 6. Heatmap plotting for a single scale
create_ggplot_heatmap <- function(results_df, title, scale_name) {
  disease_order <- scale_diseases[[scale_name]]
  diseases_in_data <- unique(results_df$Disease)
  disease_order <- disease_order[disease_order %in% diseases_in_data]
  
  results_df$Disease <- factor(results_df$Disease, levels = rev(disease_order))
  results_df$BrainRegionSimple <- factor(results_df$BrainRegionSimple, levels = brain_regions_simple)
  
  p <- ggplot(results_df, aes(x = BrainRegionSimple, y = Disease, fill = r_value)) +
    geom_tile(color = "black", linewidth = 0.5) +
    geom_text(aes(label = r_text), size = 5, color = "black", fontface = "plain", vjust = 0.3) +
    geom_text(aes(label = sig_text), size = 6, color = "black", fontface = "plain", vjust = 1.6) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, na.value = "grey90",
      limits = c(-1, 1),
      name = "Partial Correlation (ρ)"
    ) +
    labs(
      title = title,
      x = "Brain Regions",
      y = "Disease Groups",
      caption = paste("*FDR q < 0.05, **FDR q < 0.01, ***FDR q < 0.001\n",
                      "Scale:", scale_name,
                      "\nPartial correlation controlling for Age and Sex",
                      "\nCorrection: FDR corrected for 14 tests per disease per scale")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40"),
      panel.grid = element_blank()
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  
  return(p)
}

# 7. Create heatmaps for each scale
cat("\nGenerating heatmaps...\n")
heatmaps <- list()
for (scale in all_scales) {
  scale_title <- switch(scale,
                        "MOCA" = "Partial Correlation: MOCA Score with Brain Region Centile Scores",
                        "EDSS" = "Partial Correlation: EDSS Score with Brain Region Centile Scores",
                        "CASE" = "Partial Correlation: CASE Score with Brain Region Centile Scores",
                        "SDMT" = "Partial Correlation: SDMT Score with Brain Region Centile Scores",
                        "CVLT" = "Partial Correlation: CVLT Score with Brain Region Centile Scores",
                        "BVMT" = "Partial Correlation: BVMT Score with Brain Region Centile Scores",
                        "PASAT" = "Partial Correlation: PASAT Score with Brain Region Centile Scores",
                        paste("Partial Correlation:", scale, "Score with Brain Region Centile Scores"))
  if (nrow(results_by_scale[[scale]]) > 0) {
    heatmaps[[scale]] <- create_ggplot_heatmap(results_by_scale[[scale]], scale_title, scale)
  } else {
    cat(paste("Warning:", scale, "has no valid data\n"))
  }
}

# 8. Save individual heatmaps
save_heatmap <- function(heatmap_plot, filename_base) {
  width_inches <- 14
  height_inches <- 8
  png_path <- file.path("figures", paste0(filename_base, ".png"))
  tiff_path <- file.path("figures", paste0(filename_base, ".tiff"))
  pdf_path <- file.path("figures", paste0(filename_base, ".pdf"))
  
  png(png_path, width = width_inches * 300, height = height_inches * 300, res = 300)
  print(heatmap_plot)
  dev.off()
  cat("PNG saved to:", png_path, "\n")
  
  tiff(tiff_path, width = width_inches * 300, height = height_inches * 300, 
       res = 300, compression = "lzw")
  print(heatmap_plot)
  dev.off()
  cat("TIFF saved to:", tiff_path, "\n")
  
  pdf(pdf_path, width = width_inches, height = height_inches)
  print(heatmap_plot)
  dev.off()
  cat("PDF saved to:", pdf_path, "\n")
}

cat("\nSaving heatmap files...\n")
for (scale in all_scales) {
  if (!is.null(heatmaps[[scale]])) {
    save_heatmap(heatmaps[[scale]], paste0(scale, "_Partial_Correlation_Heatmap_DiseaseScale_Corrected"))
  }
}

# 9. Save the combined heatmap
save_combined_heatmap <- function(heatmap_plot, filename_base, width_inches = 14, height_inches = 24) {
  png_path <- file.path("figures", paste0(filename_base, ".png"))
  tiff_path <- file.path("figures", paste0(filename_base, ".tiff"))
  pdf_path <- file.path("figures", paste0(filename_base, ".pdf"))
  
  png(png_path, width = width_inches * 300, height = height_inches * 300, res = 300)
  print(heatmap_plot)
  dev.off()
  cat("PNG saved to:", png_path, "\n")
  
  tiff(tiff_path, width = width_inches * 300, height = height_inches * 300, 
       res = 300, compression = "lzw")
  print(heatmap_plot)
  dev.off()
  cat("TIFF saved to:", tiff_path, "\n")
  
  pdf(pdf_path, width = width_inches, height = height_inches)
  print(heatmap_plot)
  dev.off()
  cat("PDF saved to:", pdf_path, "\n")
}


cat("\nCreating a combined heatmap for all significant results...\n")

# Combine disease-scale combinations with at least one significant result
combined_data <- data.frame()
for (scale in all_scales) {
  if (nrow(results_by_scale[[scale]]) > 0) {
    # Identify diseases with at least one significant brain-region association after FDR correction
    sig_diseases <- results_by_scale[[scale]] %>%
      group_by(Disease) %>%
      filter(any(p_adj < 0.05, na.rm = TRUE)) %>%
      pull(Disease) %>%
      unique()
    
    if (length(sig_diseases) > 0) {
      scale_data <- results_by_scale[[scale]] %>%
        filter(Disease %in% sig_diseases) %>%
        mutate(Disease_Scale = paste0(Disease, " (", Scale, ")"))
      
      combined_data <- bind_rows(combined_data, scale_data)
    }
  }
}

if (nrow(combined_data) == 0) {
  cat("No significant results were found; the combined heatmap will not be generated.\n")
} else {
  y_order <- c()
  for (scale in all_scales) {
    diseases_for_scale <- scale_diseases[[scale]]
    for (disease in diseases_for_scale) {
      combo <- paste0(disease, " (", scale, ")")
      if (combo %in% combined_data$Disease_Scale) {
        y_order <- c(y_order, combo)
      }
    }
  }
  
  combined_data$Disease_Scale <- factor(combined_data$Disease_Scale, levels = rev(y_order))
  combined_data$BrainRegionSimple <- factor(combined_data$BrainRegionSimple, levels = brain_regions_simple)
  
  p_combined <- ggplot(combined_data, aes(x = BrainRegionSimple, y = Disease_Scale, fill = r_value)) +
    geom_tile(color = "black", linewidth = 0.2) +
    geom_text(aes(label = r_text), size = 4, color = "black", fontface = "plain", vjust = 0.3) +
    geom_text(aes(label = sig_text), size = 5, color = "black", fontface = "plain", vjust = 1.6) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, na.value = "grey90",
      limits = c(-1, 1),
      name = "Partial Correlation (r)"
    ) +
    labs(
      title = "Significant Partial Correlations Across All Scales",
      x = "Brain Regions",
      y = "Disease Group (Scale)",
      caption = "*FDR q < 0.05, **FDR q < 0.01, ***FDR q < 0.001\nPartial correlation controlling for Age and Sex\nCorrection: FDR corrected for 14 tests per disease per scale"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40"),
      panel.grid = element_blank()
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  

  row_height_inches <- 0.6
  n_rows <- length(unique(combined_data$Disease_Scale))
  total_height <- n_rows * row_height_inches
  

  save_combined_heatmap(
    p_combined,
    "Combined_All_Scales_Significant_Only_Merged",
    height_inches = total_height
  )

}

# 10. Save results to Excel
cat("\nSaving results to Excel file...\n")
wb <- createWorkbook()

addWorksheet(wb, "Sample_Counts")
writeData(wb, "Sample_Counts", sample_counts_all, startCol = 1, startRow = 1)
writeData(wb, "Sample_Counts", 
          data.frame(Note = "Sample counts: n = total sample size; n_valid = number of valid samples for the scale with available age and sex information."), 
          startCol = 1, startRow = nrow(sample_counts_all) + 3)

# Description of the correction method and scale-specific disease groups
correction_method <- data.frame(
  Description = c(
    "Method: partial correlation analysis controlling for age and sex.",
    "Multiple-comparison correction: disease-scale-specific FDR correction.",
    "1. For each disease and each scale, partial correlations were computed across 14 brain regions.",
    "2. Covariates: age and sex (Sex_numeric: 1 = male, 0 = female).",
    "3. FDR correction was performed separately for the 14 brain-region tests within each disease-scale combination.",
    "4. Disease groups included for each scale:",
    "   - MOCA: MS, AQP4Pos_NMOSD, MOGAD, NMDA, LGI1, GAD65",
    "   - EDSS: MS, AQP4Pos_NMOSD, MOGAD",
    "   - CASE: NMDA, LGI1, GAD65",
    "   - SDMT: MS, AQP4Pos_NMOSD",
    "   - CVLT: MS, AQP4Pos_NMOSD",
    "   - BVMT: MS, AQP4Pos_NMOSD",
    "   - PASAT: MS, AQP4Pos_NMOSD",
    "5. Each FDR correction included 14 brain-region tests.",
    "6. Significance threshold: FDR q < 0.05."
  )
)
addWorksheet(wb, "Correction_Method")
writeData(wb, "Correction_Method", correction_method, startCol = 1, startRow = 1)

for (scale in all_scales) {
  if (nrow(results_by_scale[[scale]]) > 0) {
    df_temp <- results_by_scale[[scale]] %>%
      mutate(
        r_value = round(r_value, 3),
        p_value_raw = p_value,
        p_adj_raw = p_adj,
        significance_note = case_when(
          p_adj_raw < 0.001 ~ "*** (FDR q < 0.001)",
          p_adj_raw < 0.01  ~ "** (FDR q < 0.01)",
          p_adj_raw < 0.05  ~ "* (FDR q < 0.05)",
          is.na(p_adj_raw)  ~ "NA",
          TRUE              ~ "Not significant"
        ),
        p_value = format.pval(p_value_raw, digits = 3, eps = 0.001),
        p_adj = format.pval(p_adj_raw, digits = 3, eps = 0.001)
      )
    
    scale_formatted <- df_temp %>%
      dplyr::select(dplyr::any_of(c(
        "Disease", "BrainRegionSimple", "Scale", 
        "n", "r_value", "p_value", "p_adj", "significance_note"
      )))
    
    worksheet_name <- paste0(scale, "_Partial_Results")
    addWorksheet(wb, worksheet_name)
    writeData(wb, worksheet_name, 
              data.frame(Title = paste0(scale, " partial correlation analysis controlling for age and sex with FDR correction")), 
              startCol = 1, startRow = 1)
    writeData(wb, worksheet_name, 
              data.frame(Note = "Method: partial correlation analysis with FDR correction across 14 brain-region tests."), 
              startCol = 1, startRow = 2)
    writeData(wb, worksheet_name, scale_formatted, startCol = 1, startRow = 4)
  }
}

excel_path <- file.path("results", "Partial_Correlation_Results_AllScales_DiseaseScale_Corrected.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)
cat("\nAnalysis completed.\n")
cat("Heatmap files:", file.path(output_path, "figures"), "\n")
cat("Excel results:", excel_path, "\n")
