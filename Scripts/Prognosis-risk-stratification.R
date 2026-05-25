packages_needed <- c("tidyverse", "survival", "survminer", "openxlsx", 
                     "gridExtra", "ggplot2", "patchwork", "glmnet", "caret")
for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Set paths and load data
setwd("D:/output_folder")
data_raw <- read.csv("D:/input_folder/data.csv", stringsAsFactors = FALSE)

# Define list of brain regions
brain_regions <- c(
  "Quant_lh_Hippocampus", "Quant_rh_Hippocampus",
  "Quant_lh_Amygdala", "Quant_rh_Amygdala",
  "Quant_lh_Accumbens", "Quant_rh_Accumbens",
  "Quant_lh_Anteriorcingulate", "Quant_rh_Anteriorcingulate",
  "Quant_lh_Posteriorcingulate", "Quant_rh_Posteriorcingulate",
  "Quant_lh_Entorhinal", "Quant_rh_Entorhinal",
  "Quant_lh_Orbitofrontal", "Quant_rh_Orbitofrontal"
)

# Define disease-outcome combinations to analyze
analysis_combinations <- list(
  c("MS", "Followup_relapse"),
  c("AQP4Pos_NMOSD", "Followup_relapse"),
  c("MOGAD", "Followup_relapse"),
  c("MS", "EDSS_Progress"),
  c("AQP4Pos_NMOSD", "EDSS_Progress"),
  c("MOGAD", "EDSS_Progress"),
  c("MS", "SPMS_conversion"),
  c("MS_NMOSD_Combined", "EDSS_Progress"), 
  c("Four_Disease_Combined", "Followup_relapse") 
)

# ======================== Data cleaning function ========================
clean_and_prepare_data <- function(data, diagnosis, outcome) {
  cat(paste0("\n=== Data preparation: ", diagnosis, " - ", outcome, " ===\n"))
  
  if (diagnosis == "MS_NMOSD_Combined") {
    data_subset <- data %>%
      filter(Diagnosis %in% c("MS", "AQP4Pos_NMOSD")) %>%
      mutate(
        time = as.numeric(as.character(Followup_time_month))
      )
  } else if (diagnosis == "Four_Disease_Combined") {
    # Extract combined group of multiple diseases
    data_subset <- data %>%
      filter(Diagnosis %in% c("MS", "AQP4Pos_NMOSD", "MOGAD", "NMDA", "LGI1", "GAD65")) %>%
      mutate(
        time = as.numeric(as.character(Followup_time_month))
      )
  } else {
    data_subset <- data %>%
      filter(Diagnosis == diagnosis) %>%
      mutate(
        time = as.numeric(as.character(Followup_time_month))
      )
  }
  
  # Process outcome variable
  if (outcome %in% names(data_subset)) {
    cat("Processing outcome variable '", outcome, "':\n", sep = "")
    outcome_raw <- as.numeric(as.character(data_subset[[outcome]]))
    
    status_converted <- case_when(
      is.na(outcome_raw) ~ NA_real_,
      outcome_raw == 0 ~ 0,
      outcome_raw > 0 ~ 1,
      TRUE ~ NA_real_
    )
    
    data_subset$status <- status_converted
  } else {
    stop("Error: Column '", outcome, "' does not exist")
  }
  
  # Extract brain region data and convert to numeric
  brain_data <- data_subset[, brain_regions, drop = FALSE]
  for (col in brain_regions) {
    brain_data[[col]] <- as.numeric(as.character(brain_data[[col]]))
  }
  
  # Merge survival and brain region data
  all_data <- cbind(
    data_subset[, c("time", "status")],
    brain_data
  )
  
  # Remove missing values
  complete_idx <- complete.cases(all_data)
  all_data_complete <- all_data[complete_idx, ]
  
  n_original <- nrow(all_data)
  n_complete <- nrow(all_data_complete)
  n_removed <- n_original - n_complete
  
  survival_data <- all_data_complete[, c("time", "status")]
  brain_data_clean <- all_data_complete[, brain_regions, drop = FALSE]
  
  n_total <- n_complete
  n_events <- sum(survival_data$status == 1, na.rm = TRUE)
  event_rate <- n_events / n_total * 100
  
  cat(paste0("Final sample size: ", n_total, ", Number of events: ", n_events, " (", round(event_rate, 1), "%)\n"))
  
  return(list(
    survival_data = survival_data,
    brain_data = brain_data_clean,
    n_total = n_total,
    n_events = n_events
  ))
}

# ======================== Elastic Net Cox fitting function ====================
fit_elasticnet_cox <- function(survival_data, brain_data, alpha = 0.5, s = "lambda.min") {
  x <- as.matrix(brain_data)
  y <- survival::Surv(survival_data$time, survival_data$status)
  
  set.seed(123)
  cv_fit <- cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = 10)
  
  coef_mat <- coef(cv_fit, s = s)
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  
  selected_vars <- names(coef_vec)[abs(coef_vec) > 1e-6] 
  
  return(list(
    cv_fit = cv_fit,
    coef = coef_vec,
    selected_vars = selected_vars,
    success = length(selected_vars) > 0
  ))
}

# ======================== Fit standard Cox model using selected variables ========================
fit_standard_cox <- function(survival_data, brain_data, selected_vars) {
  if (length(selected_vars) == 0) {
    return(list(model = NULL, success = FALSE))
  }
  cox_data <- cbind(survival_data, brain_data[, selected_vars, drop = FALSE])
  formula <- as.formula(paste("Surv(time, status) ~", paste(selected_vars, collapse = " + ")))
  
  tryCatch({
    model <- coxph(formula, data = cox_data)
    return(list(model = model, success = TRUE))
  }, error = function(e) {
    return(list(model = NULL, success = FALSE))
  })
}

# ============================ 10-fold nested cross-validation ============================
cross_validate_elasticnet <- function(survival_data, brain_data, alpha = 0.5, n_folds = 10) {
  set.seed(123) 
  n <- nrow(survival_data)
  
  folds <- createFolds(y = survival_data$status, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  oob_risk_scores <- rep(NA, n)
  c_indices <- c()
  selected_vars_list <- list()
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(1:n, test_idx)
    
    train_surv <- survival_data[train_idx, , drop = FALSE]
    train_brain <- brain_data[train_idx, , drop = FALSE]
    test_surv <- survival_data[test_idx, , drop = FALSE]
    test_brain <- brain_data[test_idx, , drop = FALSE]
    
    x_train <- as.matrix(train_brain)
    y_train <- Surv(train_surv$time, train_surv$status)
    
    cv_fit <- tryCatch({
      cv.glmnet(x_train, y_train, family = "cox", alpha = alpha, nfolds = min(10, sum(train_surv$status)))
    }, error = function(e) NULL)
    
    if (is.null(cv_fit)) {
      oob_risk_scores[test_idx] <- 0  
      c_indices[i] <- NA
      selected_vars_list[[i]] <- character(0)
      next
    }
    
    coef_mat <- coef(cv_fit, s = "lambda.min")
    coef_vec <- as.numeric(coef_mat)
    names(coef_vec) <- rownames(coef_mat)
    selected <- names(coef_vec)[abs(coef_vec) > 1e-6]
    selected_vars_list[[i]] <- selected
    
    if (length(selected) == 0) {
      oob_risk_scores[test_idx] <- 0
      c_indices[i] <- NA
      next
    }
    
    fit_res <- fit_standard_cox(train_surv, train_brain, selected)
    if (!fit_res$success) {
      oob_risk_scores[test_idx] <- 0
      c_indices[i] <- NA
      next
    }
    
    test_data <- cbind(test_surv, test_brain[, selected, drop = FALSE])
    lp <- predict(fit_res$model, newdata = test_data, type = "lp")
    
    oob_risk_scores[test_idx] <- lp
    
    c_index <- tryCatch({
      concordance(Surv(time, status) ~ lp, data = test_data, reverse = TRUE)$concordance
    }, error = function(e) NA)
    
    c_indices[i] <- c_index
  }
  
  mean_c <- mean(c_indices, na.rm = TRUE)
  sd_c <- sd(c_indices, na.rm = TRUE)
  cat(sprintf("  Cross-validation mean C-index = %.3f (SD = %.3f)\n", mean_c, sd_c))
  
  return(list(
    oob_risk_scores = oob_risk_scores,
    c_indices = c_indices,
    mean_c = mean_c,
    sd_c = sd_c,
    selected_vars_list = selected_vars_list
  ))
}

# ============================= Risk stratification function =============================
stratify_risk <- function(cox_data, risk_score_col = "risk_score") {
  # Uniformly use median as cutoff to avoid P-hacking suspicion
  cutoff <- median(cox_data[[risk_score_col]], na.rm = TRUE)
  cox_data$risk_group <- factor(
    ifelse(cox_data[[risk_score_col]] > cutoff, "High Risk", "Low Risk"),
    levels = c("Low Risk", "High Risk")
  )
  return(list(
    cox_data = cox_data,
    cutoff = cutoff
  ))
}

# ======================== Kaplan-Meier analysis function ========================
perform_km_analysis <- function(cox_data, diagnosis, outcome) {
  km_fit <- survfit(Surv(time, status) ~ risk_group, data = cox_data)
  logrank_test <- survdiff(Surv(time, status) ~ risk_group, data = cox_data)
  p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  p_text <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_value)))
  
  km_summary <- summary(km_fit)
  median_time_high <- NA
  median_time_low <- NA
  if (is.matrix(km_summary$table)) {
    high_risk_row <- which(grepl("High Risk", rownames(km_summary$table)))
    low_risk_row <- which(grepl("Low Risk", rownames(km_summary$table)))
    if (length(high_risk_row) > 0) median_time_high <- km_summary$table[high_risk_row[1], "median"]
    if (length(low_risk_row) > 0) median_time_low <- km_summary$table[low_risk_row[1], "median"]
  }
  
  custom_theme <- theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 22),
      axis.line = element_line(linewidth = 2.0),
      axis.ticks = element_line(linewidth = 1.5),
      legend.text = element_text(size = 22, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  km_plot <- ggsurvplot(
    km_fit,
    data = cox_data,
    pval = TRUE,
    pval.coord = c(0.1, 0.1),
    pval.size = 7,
    conf.int = TRUE,
    conf.int.alpha = 0.2,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    ncensor.plot = FALSE,
    legend.labs = c("Low Risk", "High Risk"),
    legend.title = "",
    palette = c("#006400", "#FF6F00"),
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    title = paste0(diagnosis, " - ", outcome),
    font.title = c(26, "bold"),
    font.x = c(24, "bold"),
    font.y = c(24, "bold"),
    font.tickslab = c(22, "plain"),
    font.legend = c(22, "bold"),
    size = 3.0,
    censor.size = 5,
    ggtheme = custom_theme,
    tables.theme = theme_classic() +
      theme(
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 24),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.2)
      )
  )
  
  if (!is.na(median_time_low)) {
    km_plot$plot <- km_plot$plot +
      geom_segment(aes(x = 0, y = 0.5, xend = median_time_low, yend = 0.5),
                   linetype = "dashed", color = "gray50", linewidth = 1.2) +
      geom_segment(aes(x = median_time_low, y = 0.5, xend = median_time_low, yend = 0),
                   linetype = "dashed", color = "gray50", linewidth = 1.2)
  }
  if (!is.na(median_time_high)) {
    km_plot$plot <- km_plot$plot +
      geom_segment(aes(x = 0, y = 0.5, xend = median_time_high, yend = 0.5),
                   linetype = "dashed", color = "gray50", linewidth = 1.2) +
      geom_segment(aes(x = median_time_high, y = 0.5, xend = median_time_high, yend = 0),
                   linetype = "dashed", color = "gray50", linewidth = 1.2)
  }
  
  km_plot$plot <- km_plot$plot +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))
  
  return(list(
    km_fit = km_fit,
    km_plot = km_plot,
    p_value = p_value,
    p_text = p_text,
    median_time_high = median_time_high,
    median_time_low = median_time_low
  ))
}

# ======================== Main analysis function ============================
analyze_combination <- function(data, diagnosis, outcome, 
                                alpha = 0.5,
                                output_dir = "survival_analysis_results_elasticnet") {
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Starting analysis:", diagnosis, " - ", outcome, "\n")
  cat(strrep("=", 60), "\n")
  
  # Create output directories
  dirs_to_create <- c(
    output_dir,
    file.path(output_dir, "KM_plots"),
    file.path(output_dir, "ElasticNet_results"),
    file.path(output_dir, "CV_results"),
    file.path(output_dir, "Summary_results"),
    file.path(output_dir, "Selected_Variables")
  )
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 1. Data preparation
  prepared_data <- clean_and_prepare_data(data, diagnosis, outcome)
  if (prepared_data$n_total < 20) {
    cat("Warning: Sample size too small, skipping analysis\n")
    return(NULL)
  }
  if (prepared_data$n_events < 10) {
    cat("Warning: Few events, interpret results with caution\n")
  }
  
  # 2. 10-fold nested CV prediction (to obtain OOB predictions without leakage for KM)
  cat("\n--- Perform 10-fold nested CV prediction (eliminate data leakage) ---\n")
  cv_results <- cross_validate_elasticnet(
    prepared_data$survival_data, 
    prepared_data$brain_data,
    alpha = alpha,
    n_folds = 10
  )
  
  cv_df <- data.frame(
    Fold = 1:10,
    C_index = cv_results$c_indices,
    Selected_Vars_Count = sapply(cv_results$selected_vars_list, length)
  )
  write.csv(cv_df, file.path(output_dir, "CV_results", paste0("CV_", diagnosis, "_", outcome, ".csv")), row.names = FALSE)
  
  # 3. Fit final Elastic Net model on full dataset
  cat("\n--- Fit final Elastic Net model on full dataset ---\n")
  enet_res <- fit_elasticnet_cox(prepared_data$survival_data, prepared_data$brain_data, alpha = alpha)
  
  final_cindex <- NA
  num_selected_vars <- 0
  
  if (enet_res$success) {
    selected_vars <- enet_res$selected_vars
    num_selected_vars <- length(selected_vars)
    cat("Variables selected by full data model (", num_selected_vars, "): ", paste(selected_vars, collapse = ", "), "\n", sep = "")
    
    write.csv(data.frame(Variable = selected_vars, Coefficient = enet_res$coef[selected_vars]), 
              file.path(output_dir, "Selected_Variables", paste0("SelectedVars_", diagnosis, "_", outcome, ".csv")), row.names = FALSE)
    
    final_fit <- fit_standard_cox(prepared_data$survival_data, prepared_data$brain_data, selected_vars)
    if (final_fit$success) {
      cox_summary <- summary(final_fit$model)
      final_cindex <- cox_summary$concordance[1]
      
      # Save full data Cox coefficients
      cox_df <- as.data.frame(cox_summary$coefficients)
      cox_df$Variable <- rownames(cox_df)
      cox_df <- cox_df[, c("Variable", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")]
      names(cox_df) <- c("Variable", "Coefficient", "HR", "SE", "z", "P_value")
      
      ci_df <- as.data.frame(cox_summary$conf.int)
      cox_df$HR_95CI_Lower <- ci_df$`lower .95`
      cox_df$HR_95CI_Upper <- ci_df$`upper .95`
      cox_df$HR_95CI <- sprintf("%.2f (%.2f-%.2f)", cox_df$HR, cox_df$HR_95CI_Lower, cox_df$HR_95CI_Upper)
      cox_df$Significance <- ifelse(cox_df$P_value < 0.001, "***",
                                    ifelse(cox_df$P_value < 0.01, "**",
                                           ifelse(cox_df$P_value < 0.05, "*", "NS")))
      write.csv(cox_df, file.path(output_dir, "ElasticNet_results", paste0("Cox_Detailed_", diagnosis, "_", outcome, ".csv")), row.names = FALSE)
    }
  } else {
    cat("Elastic Net failed to select any variables on the full dataset. Still generating KM plot based on CV OOB scores.\n")
  }
  
  # 4. Assemble and clean data for KM: use non-leaking predicted values (OOB) from CV
  cox_data <- prepared_data$survival_data
  cox_data$risk_score <- cv_results$oob_risk_scores
  
  valid_idx <- !is.na(cox_data$risk_score)
  cox_data <- cox_data[valid_idx, ]
  
  # 5. Uniformly use median for risk stratification
  stratified <- stratify_risk(cox_data, risk_score_col = "risk_score")
  cox_data_risk <- stratified$cox_data
  
  # 6. Generate KM curves and test
  km_results <- perform_km_analysis(cox_data_risk, diagnosis, outcome)
  
  # --- Calculate overall hazard ratio (Group HR) between high and low risk groups ---
  group_cox <- coxph(Surv(time, status) ~ risk_group, data = cox_data_risk)
  group_sum <- summary(group_cox)
  group_hr <- group_sum$conf.int[1, "exp(coef)"]
  group_lower <- group_sum$conf.int[1, "lower .95"]
  group_upper <- group_sum$conf.int[1, "upper .95"]
  group_p <- group_sum$coefficients[1, "Pr(>|z|)"]
  # -----------------------------------------------------------
  
  ggsave(file.path(output_dir, "KM_plots", paste0("KM_", diagnosis, "_", outcome, ".png")), 
         print(km_results$km_plot), width = 10, height = 10, dpi = 300, bg = "white")
  
  # 7. Summarize and return
  summary_df <- data.frame(
    Diagnosis = diagnosis,
    Outcome = outcome,
    N_Total = prepared_data$n_total,
    N_Events = prepared_data$n_events,
    Event_Rate = round(prepared_data$n_events / prepared_data$n_total * 100, 1),
    Alpha = alpha,
    N_Selected_Vars = num_selected_vars,
    CV_Mean_Cindex = round(cv_results$mean_c, 3),
    CV_SD_Cindex = round(cv_results$sd_c, 3),
    Final_Cindex = round(final_cindex, 3),
    Group_HR = round(group_hr, 2),                             
    Group_HR_95CI = sprintf("%.2f-%.2f", group_lower, group_upper), 
    Group_Cox_P = group_p,                                     
    KM_P_Value = km_results$p_value,
    Median_Time_Low_Risk = round(km_results$median_time_low, 1),
    Median_Time_High_Risk = round(km_results$median_time_high, 1),
    Time_Difference = round(km_results$median_time_high - km_results$median_time_low, 1),
    Risk_Cutoff = round(stratified$cutoff, 3),
    KM_Significance = ifelse(km_results$p_value < 0.001, "***",
                             ifelse(km_results$p_value < 0.01, "**",
                                    ifelse(km_results$p_value < 0.05, "*", "NS")))
  )
  write.csv(summary_df, file.path(output_dir, "Summary_results", paste0("Summary_", diagnosis, "_", outcome, ".csv")), row.names = FALSE)
  
  return(list(
    summary_df = summary_df,
    cv_results = cv_results,
    km_results = km_results
  ))
}

# ======================== Batch analysis function ========================
batch_analyze <- function(data, analysis_combinations, alpha = 0.5,
                          output_dir = "survival_analysis_results_elasticnet") {
  all_results <- list()
  summary_table <- data.frame()
  
  for (i in seq_along(analysis_combinations)) {
    diagnosis <- analysis_combinations[[i]][1]
    outcome <- analysis_combinations[[i]][2]
    
    result <- analyze_combination(data, diagnosis, outcome, 
                                  alpha = alpha,
                                  output_dir = output_dir)
    
    if (!is.null(result)) {
      all_results[[paste0(diagnosis, "_", outcome)]] <- result
      summary_table <- rbind(summary_table, result$summary_df)
    } 
  }
  
  return(list(
    all_results = all_results,
    summary_table = summary_table
  ))
}

# ======================== Create Excel report ========================
create_excel_report <- function(results, output_dir = "survival_analysis_results_elasticnet") {
  cat("\nCreating Excel report...\n")
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", results$summary_table, startRow = 1, startCol = 1)
  setColWidths(wb, "Summary", cols = 1:ncol(results$summary_table), widths = "auto")
  
  significant_results <- results$summary_table %>%
    filter(!is.na(KM_P_Value), KM_P_Value < 0.05) %>%
    dplyr::select(Diagnosis, Outcome, N_Total, N_Events, Group_HR, Group_HR_95CI, Group_Cox_P, KM_P_Value, KM_Significance,
                  Median_Time_Low_Risk, Median_Time_High_Risk, Time_Difference, 
                  CV_Mean_Cindex, Final_Cindex)
  if (nrow(significant_results) > 0) {
    addWorksheet(wb, "Significant_Results")
    writeData(wb, "Significant_Results", significant_results)
  }
  
  addWorksheet(wb, "Analysis_Notes")
  notes_df <- data.frame(
    Item = c("Analysis Method", "Zero Data Leakage", "Risk Cutoff", "C-index Interpretation", "Group HR"),
    Description = c(
      "Elastic Net feature selection embedded in Nested CV + Kaplan-Meier survival analysis.",
      "The risk scores used for KM curves are strictly Out-of-bag predictions derived from independent test folds, completely preventing data leakage.",
      "To ensure robustness, patients are dichotomized into high/low risk groups using the median risk score as the standard cutoff.",
      "C-index from full-data models can be over-optimistic; rely primarily on CV_Mean_Cindex to evaluate out-of-sample generalization.",
      "The overall Hazard Ratio (High vs. Low risk group) based on OOB risk scores, reflecting the relative risk of progression/relapse."
    )
  )
  writeData(wb, "Analysis_Notes", notes_df)
  
  excel_file <- file.path(output_dir, "survival_analysis_elasticnet_results.xlsx")
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  return(excel_file)
}

# ======================== Create combined KM plot ========================
create_combined_km_plot <- function(output_dir = "survival_analysis_results_elasticnet") {
  km_files <- list.files(file.path(output_dir, "KM_plots"), 
                         pattern = "\\.png$", full.names = TRUE)
  if (length(km_files) == 0) return(NULL)
  
  if (!require(png, quietly = TRUE)) install.packages("png")
  library(png)
  
  km_images <- lapply(km_files, readPNG)
  n_plots <- length(km_images)
  n_cols <- 3
  n_rows <- ceiling(n_plots / n_cols)
  
  combined_file <- file.path(output_dir, "combined_KM_plots.png")
  png(combined_file, width = 3800, height = 2600 * n_rows / 2, res = 300)
  
  layout_matrix <- matrix(1:(n_cols * n_rows), nrow = n_rows, ncol = n_cols, byrow = TRUE)
  layout(layout_matrix)
  par(mar = c(3, 3, 4, 2))
  
  for (i in 1:n_plots) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), 
         xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    rasterImage(km_images[[i]], 0, 0, 1, 1)
    title_text <- gsub("KM_", "", gsub("\\.png$", "", basename(km_files[i])))
    title_text <- gsub("_", " - ", title_text)
    title(main = title_text, cex.main = 1.6, line = 0.5)
  }
  
  dev.off()
}

# ======================== Generate text report ========================
generate_report <- function(results, output_dir = "survival_analysis_results_elasticnet") {
  report_file <- file.path(output_dir, "analysis_report.txt")
  sink(report_file)
  cat("Limbic System Brain Region Survival Analysis Report (Elastic Net - Zero Leakage)\n")
  cat("Generated on: ", as.character(Sys.time()), "\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n\n")
  cat("This code leverages Nested CV to produce out-of-bag blind test risk scores, ensuring the generated Kaplan-Meier curves are free of data leakage.\n\n")
  
  cat("Summary of significant results (KM p < 0.05)\n")
  cat(paste0(rep("-", 60), collapse = ""), "\n")
  significant <- results$summary_table %>% filter(!is.na(KM_P_Value), KM_P_Value < 0.05) %>% arrange(KM_P_Value)
  if (nrow(significant) > 0) {
    for (i in 1:nrow(significant)) {
      row <- significant[i, ]
      cat(sprintf("%s - %s: KM p = %.4f%s, Group HR = %s, CV C-index = %.3f\n", 
                  row$Diagnosis, row$Outcome, row$KM_P_Value, row$KM_Significance, row$Group_HR_95CI, row$CV_Mean_Cindex))
    }
  } else {
    cat("No significant results found (p < 0.05)\n")
  }
  sink()
}

# ======================== Main execution ========================
cat("=== Limbic System Brain Region Survival Analysis (Elastic Net + OOB Zero Leakage) Started ===\n")

# Check required data columns
required_cols <- c("Diagnosis", "Followup_time_month", 
                   "Followup_relapse", "EDSS_Progress", "SPMS_conversion", brain_regions)
missing_cols <- setdiff(required_cols, names(data_raw))
if (length(missing_cols) > 0) stop("Missing following columns: ", paste(missing_cols, collapse = ", "))

# Batch analysis
alpha <- 0.5   
output_dir <- "survival_analysis_results_elasticnet"
results <- batch_analyze(data_raw, analysis_combinations, alpha = alpha, output_dir = output_dir)

# Generate reports
excel_file <- create_excel_report(results, output_dir)
create_combined_km_plot(output_dir)
generate_report(results, output_dir)

cat("\n=== Finished! Results written to ", output_dir, "===\n")