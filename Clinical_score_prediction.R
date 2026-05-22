library(glmnet)      
library(e1071)       
library(ggplot2)     
library(openxlsx)    


# ==================== Set paths and parameters ====================
input_file <- "D:/input_folder/data.csv"
output_dir <- "D:/output_folder"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define disease-scale combinations
combinations <- list(
  c("MS", "MOCA"),
  c("AQP4Pos_NMOSD", "MOCA"),
  c("MOGAD", "MOCA"),
  c("LGI1", "MOCA"),
  c("NMDA", "MOCA"),
  c("MS", "SDMT"),
  c("MS", "CVLT"),
  c("AQP4Pos_NMOSD", "CVLT"),
  c("MS", "BVMT"),
  c("AQP4Pos_NMOSD", "BVMT"),
  c("MS", "EDSS"),
  c("AQP4Pos_NMOSD", "EDSS"),
  c("MOGAD", "EDSS"),
  c("LGI1", "CASE"),
  c("NMDA", "CASE")
)

results_summary <- data.frame()

# ==================== Load data ====================
data_all <- read.csv(input_file, stringsAsFactors = FALSE)


# Automatically identify brain-region centile-score features starting with "Quant_"
feature_cols <- grep("^Quant_", names(data_all), value = TRUE)
cat("Found", length(feature_cols), "brain-region feature columns:\n")
print(feature_cols)

# ==================== Process each disease-scale combination ====================
for (comb in combinations) {
  diagnosis <- comb[1]
  scale <- comb[2]
  
  cat("\n===== Processing", diagnosis, "-", scale, "=====\n")
  
  # ---------- Filter data by diagnosis ----------
    subset <- data_all[data_all$Diagnosis == diagnosis, ]
  
  if (nrow(subset) == 0) {
    cat("Warning: no data for diagnosis", diagnosis, "; skipped\n")
    next
  }
  
  # Check whether the scale column exists
  if (!(scale %in% names(subset))) {
    cat("Warning: scale", scale, "does not exist in the data; skipped\n")
    next
  }
  
  subset <- subset[!is.na(subset[[scale]]), ]
  if (nrow(subset) == 0) {
    cat("Warning:", scale, "is entirely NA in diagnosis", diagnosis, "; skipped\n")
    next
  }
  
  # ---------- Process covariates: age and sex ----------
  # Check whether Age and Sex columns exist
  if (!("Age" %in% names(subset)) || !("Sex" %in% names(subset))) {
    cat("Warning: Age or Sex column is missing; skipped\n")
    next
  }
  # Remove samples with missing Age or Sex values
  subset <- subset[!is.na(subset$Age) & !is.na(subset$Sex), ]
  if (nrow(subset) == 0) {
    cat("Warning: Age or Sex contains missing values; skipped\n")
    next
  }
  # Convert Sex to numeric: Male = 1, Female = 0
  subset$Sex_num <- ifelse(subset$Sex == "Male", 1, 0)
  
  # Construct the feature matrix including brain-region features and covariates
  all_features <- c(feature_cols, "Age", "Sex_num")
  X <- as.matrix(subset[, all_features])  
  y <- subset[[scale]]
  feature_names <- colnames(X) 
  
  n_samples <- nrow(X)
  cat("Valid sample size:", n_samples, "\n")
  
  if (n_samples < 5) {
    cat("Sample size is too small for reliable prediction; skipped\n")
    next
  }
  
  # Determine the number of cross-validation folds
  if (n_samples < 10) {
    n_folds <- n_samples
    shuffle <- FALSE
    cat("Sample size < 10; using leave-one-out cross-validation (n_folds =", n_folds, ")\n")
  } else {
    n_folds <- 10
    shuffle <- TRUE
  }
  
  # Generate cross-validation folds
  set.seed(42)
  if (shuffle) {
    fold_id <- sample(rep(1:n_folds, length.out = n_samples))
  } else {
    fold_id <- 1:n_samples 
  }
  
  y_true_all <- c()
  y_pred_all <- c()
  
  # Cross-validation loop
  for (fold in 1:n_folds) {
    train_idx <- which(fold_id != fold)
    test_idx  <- which(fold_id == fold)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]
    y_test  <- y[test_idx]
    
    # 1. Standardize features using the training set
    mean_train <- apply(X_train, 2, mean)
    sd_train   <- apply(X_train, 2, sd)
    sd_train[sd_train == 0] <- 1
    X_train_scaled <- scale(X_train, center = mean_train, scale = sd_train)
    X_test_scaled  <- scale(X_test, center = mean_train, scale = sd_train)
    
    # 2. Elastic net feature selection (alpha = 0.5)
    penalty_factors <- rep(1, ncol(X_train_scaled))
    age_idx <- which(colnames(X_train_scaled) == "Age")
    sex_idx <- which(colnames(X_train_scaled) == "Sex_num")
    penalty_factors[c(age_idx, sex_idx)] <- 0 
    
    inner_nfolds <- min(5, nrow(X_train_scaled)) 
    
    cv_enet <- cv.glmnet(X_train_scaled, y_train, alpha = 0.5,
                         nfolds = inner_nfolds, standardize = FALSE,
                         penalty.factor = penalty_factors)
    best_lambda <- cv_enet$lambda.min
    enet_model <- glmnet(X_train_scaled, y_train, alpha = 0.5, lambda = best_lambda,
                         standardize = FALSE, penalty.factor = penalty_factors)
    
    coef_enet <- as.matrix(coef(enet_model))
    selected_features <- which(coef_enet[-1, 1] != 0)
    
    
    if (length(selected_features) == 0) {
      selected_features <- 1:ncol(X_train_scaled)
      cat("Warning: no features were selected; using all features\n")
    }
    
    # 3. Train SVR using selected features
    X_train_selected <- X_train_scaled[, selected_features, drop = FALSE]
    X_test_selected  <- X_test_scaled[, selected_features, drop = FALSE]
    
    svr_model <- svm(X_train_selected, y_train, type = "eps-regression", kernel = "radial")
    
    # 4. Predict
    y_pred <- predict(svr_model, X_test_selected)
    
    y_true_all <- c(y_true_all, y_test)
    y_pred_all <- c(y_pred_all, y_pred)
  }
  
  # ============= Compute performance metrics and confidence intervals ================
  # Run Pearson correlation test and extract the confidence interval
  cor_test_res <- cor.test(y_true_all, y_pred_all, method = "pearson")
  r <- cor_test_res$estimate
  p_value <- cor_test_res$p.value
  
  # Extract the 95% confidence interval
  ci_lower <- cor_test_res$conf.int[1]
  ci_upper <- cor_test_res$conf.int[2]
  r_95ci <- sprintf("[%.3f, %.3f]", ci_lower, ci_upper) 
  
  r2 <- 1 - sum((y_true_all - y_pred_all)^2) / sum((y_true_all - mean(y_true_all))^2)
  rmse <- sqrt(mean((y_true_all - y_pred_all)^2))
  mae <- mean(abs(y_true_all - y_pred_all))
  n <- length(y_true_all)
  
  cat("r =", round(r, 3), "95%CI:", r_95ci, ", p =", format(p_value, scientific = TRUE),
      ", R2 =", round(r2, 3), ", RMSE =", round(rmse, 3),
      ", MAE =", round(mae, 3), ", N =", n, "\n")
  

  mean_all <- apply(X, 2, mean)
  sd_all <- apply(X, 2, sd)
  sd_all[sd_all == 0] <- 1
  X_scaled_all <- scale(X, center = mean_all, scale = sd_all)
  
  penalty_factors_all <- rep(1, ncol(X_scaled_all))
  age_idx_all <- which(colnames(X_scaled_all) == "Age")
  sex_idx_all <- which(colnames(X_scaled_all) == "Sex_num")
  penalty_factors_all[c(age_idx_all, sex_idx_all)] <- 0
  
  cv_all_nfolds <- min(5, nrow(X_scaled_all))
  
  cv_all <- cv.glmnet(X_scaled_all, y, alpha = 0.5, nfolds = cv_all_nfolds, standardize = FALSE,
                      penalty.factor = penalty_factors_all)
  best_lambda_all <- cv_all$lambda.min
  final_model <- glmnet(X_scaled_all, y, alpha = 0.5, lambda = best_lambda_all, standardize = FALSE,
                        penalty.factor = penalty_factors_all)
  
  coef_final <- as.matrix(coef(final_model))
  selected_indices <- which(coef_final[-1, 1] != 0)
  
  
  selected_feature_names <- feature_names[selected_indices]
  selected_features_str <- paste(selected_feature_names, collapse = ", ")
  # =========================================================
  
  # ==================== Record results ===========================
  results_summary <- rbind(results_summary, data.frame(
    Diagnosis = diagnosis,
    Scale = scale,
    Pearson_r = r,
    Pearson_r_95CI = r_95ci, 
    p_value = p_value,
    R2 = r2,
    RMSE = rmse,
    MAE = mae,
    N = n,
    Selected_Features = selected_features_str,
    stringsAsFactors = FALSE
  ))
  
  # ==================== Plotting: set point color by diagnosis ====================
  if (diagnosis == "MS") {
    point_color <- "#0072B5" 
  } else if (diagnosis == "AQP4Pos_NMOSD") {
    point_color <- "#E64B35"  
  } else if (diagnosis == "MOGAD") {
    point_color <- "#228B22"   
  } else if (diagnosis %in% c("LGI1", "NMDA")) {
    point_color <- "#6A0DAD"   
  } else {
    point_color <- "#3C5488"  
  }
  
  plot_data <- data.frame(Observed = y_true_all, Predicted = y_pred_all)
  
  min_val <- min(c(plot_data$Observed, plot_data$Predicted), na.rm = TRUE)
  max_val <- max(c(plot_data$Observed, plot_data$Predicted), na.rm = TRUE)
  padding <- (max_val - min_val) * 0.02 
  axis_min <- min_val - padding
  axis_max <- max_val + padding
  
  p_base <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                linewidth = 1, color = "gray50") +
    geom_point(size = 4.5, alpha = 0.85, color = point_color) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "red", fill = "gray70", alpha = 0.3) + 
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max), expand = FALSE) +
    labs(title = paste0(diagnosis, " - ", scale),
         x = paste0("Observed ", scale),
         y = paste0("Predicted ", scale)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18), 
      panel.grid = element_blank(),
      
      axis.text = element_text(size = 18, color = "black", face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      panel.border = element_rect(linewidth = 2.5, color = "black"),
      axis.ticks = element_line(linewidth = 2, color = "black"),
      axis.ticks.length = unit(0.25, "cm"), 
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15) 
    )
  
  label_text <- paste0(
    "n = ", n,
    "\nr = ", sprintf("%.3f", r),
    "\np = ", formatC(p_value, format = "e", digits = 2),
    "\nR² = ", sprintf("%.3f", r2),
    "\nRMSE = ", sprintf("%.3f", rmse),
    "\nMAE = ", sprintf("%.3f", mae)
  )
  
  p_with_anno <- p_base +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.2,
             label = label_text, size = 5.5, fontface = "bold") 
  
  # Save annotated TIFF figure
  base_filename <- file.path(output_dir, paste0(diagnosis, "_", scale, "_prediction"))
  ggsave(paste0(base_filename, "_with_anno.tiff"), plot = p_with_anno,
         width = 7.2, height = 7.2, dpi = 300, device = "tiff", compression = "lzw")
  # =====================================================================
}

# ==================== Save result summary ====================
if (nrow(results_summary) > 0) {
  results_summary$p_fdr <- p.adjust(results_summary$p_value, method = "fdr")
  results_summary <- results_summary[order(results_summary$Diagnosis, results_summary$Scale), ]
  
  output_excel <- file.path(output_dir, "prediction_summary.xlsx")
  write.xlsx(results_summary, output_excel, rowNames = FALSE)
  cat("\nResult summary saved to:", output_excel, "\n")
} else {
  cat("\nNo results were generated. Please check the data.\n")
}

cat("\nAll tasks completed.\n")