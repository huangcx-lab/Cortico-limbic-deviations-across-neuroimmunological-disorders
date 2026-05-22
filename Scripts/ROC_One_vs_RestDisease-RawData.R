library(tidyverse)
library(glmnet)    
library(caret)      
library(e1071)      
library(pROC)     
library(MatchIt)    
library(openxlsx)  

############################################################
## 1. Load data, preprocess, and split train/test sets
############################################################
data_path <- "D:/input_folder/data.csv"

df <- read.csv(data_path)

df$Diagnosis <- as.factor(df$Diagnosis)

# Convert Sex to numeric for matching and pre-processing
df$Sex <- as.character(df$Sex)
df$Sex[df$Sex == "Male"] <- 1
df$Sex[df$Sex == "Female"] <- 0
df$Sex <- as.numeric(df$Sex)

if ("Site_ZZZ" %in% colnames(df)) {
  df$Site_ZZZ <- as.character(df$Site_ZZZ)
}

test_sites <- c("ExternalSite_01")

df_train <- df %>% filter(!Site_ZZZ %in% test_sites)
df_test  <- df %>% filter(Site_ZZZ %in% test_sites)

############################################################
## 2. Define feature set
############################################################
brain_features <- c(
  "Raw_lh_Hippocampus","Raw_rh_Hippocampus",
  "Raw_lh_Amygdala","Raw_rh_Amygdala",
  "Raw_lh_Accumbens","Raw_rh_Accumbens",
  "Raw_lh_Anteriorcingulate","Raw_rh_Anteriorcingulate",
  "Raw_lh_Posteriorcingulate","Raw_rh_Posteriorcingulate",
  "Raw_lh_Entorhinal","Raw_rh_Entorhinal",
  "Raw_lh_Orbitofrontal","Raw_rh_Orbitofrontal"
)

features_all <- brain_features

############################################################
## 3. Define matching, test-set filtering, and support-range functions
############################################################


get_binary_data <- function(data, target_group) {
  group_defs <- list(
    MS = c("MS"),
    NMOSD = c("AQP4Pos_NMOSD"),
    MOGAD = c("MOGAD"),
    AE = c("NMDA", "LGI1", "GABA", "GAD65", "CASPR2")
  )
  
  pos_diagnoses <- group_defs[[target_group]]
  neg_diagnoses <- unlist(group_defs[names(group_defs) != target_group])
  
  sub <- data %>%
    filter(Diagnosis %in% c(pos_diagnoses, neg_diagnoses)) %>%
    mutate(label = ifelse(Diagnosis %in% pos_diagnoses, 1, 0))
  
  sub <- sub[complete.cases(sub[, c("label", "Age", "Sex")]), ]
  return(sub)
}

match_data <- function(data, target_group) {
  sub <- get_binary_data(data, target_group)
  
  if (nrow(sub) == 0 || length(unique(sub$label)) < 2) {
    return(sub[0, , drop = FALSE])
  }
  
  m <- tryCatch(
    matchit(
      label ~ Age,
      data = sub,
      method = "nearest",
      exact = ~ Sex,
      caliper = c(Age = 3),
      std.caliper = FALSE,
      ratio = 1
    ),
    error = function(e) NULL
  )
  
  if (is.null(m)) {
    return(sub[0, , drop = FALSE])
  } else {
    return(match.data(m))
  }
}

filter_test_data <- function(data, target_group) {
  return(get_binary_data(data, target_group))
}

apply_train_support_range <- function(data, train_age_min, train_age_max) {
  data %>%
    filter(Age >= train_age_min & Age <= train_age_max)
}

############################################################
## 4. Core function
############################################################
run_disease_classifier <- function(d_train, d_test, disease_name, use_seed = NULL, record_seeds = FALSE) {
  
  seed_records <- list()
  if (record_seeds) {
    if (is.null(use_seed)) use_seed <- sample(1:1000000, 1)
    seed_records[["global_seed"]] <- use_seed
  }
  if (!is.null(use_seed)) set.seed(use_seed)
  
  sub_all_train <- filter_test_data(d_train, disease_name)
  sub_test_raw  <- filter_test_data(d_test, disease_name)
  
  if (nrow(sub_all_train) == 0) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("No available training samples for", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  y_train_all <- sub_all_train$label
  
  n_pos <- sum(y_train_all == 1)
  n_neg <- sum(y_train_all == 0)
  n_min <- min(n_pos, n_neg)
  
  if (n_min < 5) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Too few samples for stable CV:", disease_name, "(minority class <", 5, ")"),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  if (n_min >= 10) {
    k_cv <- 10
  } else {
    k_cv <- 5
  }
  
  if (!is.null(use_seed)) set.seed(use_seed + 999)
  folds <- createFolds(y_train_all, k = k_cv)
  
  preds_cv_int <- rep(NA, length(y_train_all))
  numeric_cols <- names(sub_all_train)[sapply(sub_all_train, is.numeric)]
  
  cv_failed <- FALSE
  cv_fail_reason <- NULL
  
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(y_train_all), test_idx)
    
    cv_train_raw <- sub_all_train[train_idx, ]
    cv_test_raw  <- sub_all_train[test_idx, ]
    
    fold_medians <- sapply(cv_train_raw[, numeric_cols, drop = FALSE], function(x) median(x, na.rm = TRUE))
    
    cv_train_imp <- cv_train_raw
    cv_test_imp  <- cv_test_raw
    
    for (col in numeric_cols) {
      if (any(is.na(cv_train_imp[[col]]))) cv_train_imp[[col]][is.na(cv_train_imp[[col]])] <- fold_medians[[col]]
      if (any(is.na(cv_test_imp[[col]])))  cv_test_imp[[col]][is.na(cv_test_imp[[col]])]  <- fold_medians[[col]]
    }
    
    cv_train_matched <- match_data(cv_train_imp, disease_name)
    
    if (nrow(cv_train_matched) < 10) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: matched train samples < 10")
      break
    }
    
    if (length(unique(cv_train_matched$label)) < 2) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: matched train labels < 2 classes")
      break
    }
    
    x_train_mat <- model.matrix(label ~ ., cv_train_matched[, c(features_all, "label")])[,-1]
    y_train_mat <- cv_train_matched$label
    
    if (record_seeds && !is.null(use_seed)) {
      lasso_seed <- use_seed + 1000 + i
      set.seed(lasso_seed)
      seed_records[[paste0("lasso_fold_", i, "_seed")]] <- lasso_seed
    }
    
    cv.lasso <- tryCatch(
      cv.glmnet(x_train_mat, y_train_mat, family = "binomial", alpha = 1, nfolds = 5),
      error = function(e) NULL
    )
    
    if (is.null(cv.lasso)) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: cv.glmnet error")
      break
    }
    
    coef_lasso <- coef(cv.lasso, s = "lambda.min")
    sel_feat <- rownames(coef_lasso)[coef_lasso[,1] != 0]
    sel_feat <- setdiff(sel_feat, "(Intercept)")
    if (length(sel_feat) == 0) sel_feat <- features_all
    
    svm_train_df <- cv_train_matched[, c("label", sel_feat), drop = FALSE]
    svm_train_df$label <- as.factor(svm_train_df$label)
    
    if (record_seeds && !is.null(use_seed)) {
      svm_seed <- use_seed + 3000 + i
      set.seed(svm_seed)
      seed_records[[paste0("svm_fold_", i, "_seed")]] <- svm_seed
    }
    
    svm_model <- tryCatch(
      svm(label ~ ., data = svm_train_df, kernel = "radial", scale = TRUE, probability = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(svm_model)) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: SVM training error")
      break
    }
    
    pred_int <- tryCatch(
      predict(svm_model, newdata = cv_test_imp[, sel_feat, drop = FALSE], probability = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(pred_int)) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: SVM prediction error")
      break
    }
    
    prob_mat <- attr(pred_int, "probabilities")
    if (is.null(prob_mat) || !("1" %in% colnames(prob_mat))) {
      cv_failed <- TRUE
      cv_fail_reason <- paste("Fold", i, "failed: probability output invalid")
      break
    }
    
    preds_cv_int[test_idx] <- prob_mat[, "1"]
  }
  
  if (cv_failed || any(is.na(preds_cv_int))) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Internal CV unstable for", disease_name, "-", cv_fail_reason),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  roc_int <- roc(y_train_all, preds_cv_int, quiet = TRUE, ci = TRUE, direction = "<")
  coords_int <- coords(
    roc_int, "best",
    ret = c("threshold", "specificity", "sensitivity", "accuracy"),
    transpose = FALSE, best.method = "youden"
  )
  if (is.data.frame(coords_int) && nrow(coords_int) > 1) coords_int <- coords_int[1, ]
  
  final_medians <- sapply(sub_all_train[, numeric_cols, drop = FALSE], function(x) median(x, na.rm = TRUE))
  final_train_imp <- sub_all_train
  
  for (col in numeric_cols) {
    if (any(is.na(final_train_imp[[col]]))) final_train_imp[[col]][is.na(final_train_imp[[col]])] <- final_medians[[col]]
  }
  
  final_train_matched <- match_data(final_train_imp, disease_name)
  n_matched <- nrow(final_train_matched)
  
  if (n_matched < 60 || length(unique(final_train_matched$label)) < 2) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Matched samples too few or one class only for final model:", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  x_final_mat <- model.matrix(label ~ ., final_train_matched[, c(features_all, "label")])[,-1]
  y_final_mat <- final_train_matched$label
  
  if (record_seeds && !is.null(use_seed)) {
    final_lasso_seed <- use_seed + 5000
    set.seed(final_lasso_seed)
    seed_records[["final_lasso_seed"]] <- final_lasso_seed
  }
  
  final_cv_lasso <- tryCatch(
    cv.glmnet(x_final_mat, y_final_mat, family = "binomial", alpha = 1, nfolds = 5),
    error = function(e) NULL
  )
  
  if (is.null(final_cv_lasso)) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Final cv.glmnet failed for", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  coef_final_lasso <- coef(final_cv_lasso, s = "lambda.min")
  final_sel_feat <- rownames(coef_final_lasso)[coef_final_lasso[,1] != 0]
  final_sel_feat <- setdiff(final_sel_feat, "(Intercept)")
  if (length(final_sel_feat) == 0) final_sel_feat <- features_all
  
  final_svm_train_df <- final_train_matched[, c("label", final_sel_feat), drop = FALSE]
  final_svm_train_df$label <- as.factor(final_svm_train_df$label)
  
  if (record_seeds && !is.null(use_seed)) {
    final_svm_seed <- use_seed + 7000
    set.seed(final_svm_seed)
    seed_records[["final_svm_seed"]] <- final_svm_seed
  }
  
  final_svm_model <- tryCatch(
    svm(label ~ ., data = final_svm_train_df, kernel = "radial", scale = TRUE, probability = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(final_svm_model)) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Final SVM training failed for", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  pred_final_int <- tryCatch(
    predict(final_svm_model, newdata = final_train_imp[, final_sel_feat, drop = FALSE], probability = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(pred_final_int)) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Final SVM internal prediction failed for thresholding:", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  prob_final_int <- attr(pred_final_int, "probabilities")
  if (is.null(prob_final_int) || !("1" %in% colnames(prob_final_int))) {
    return(list(
      int = NULL, ext = NULL,
      reason = paste("Final SVM probability output invalid for thresholding:", disease_name),
      full_preds_train = NULL, full_preds_test = NULL,
      seed_records = if(record_seeds) seed_records else NULL
    ))
  }
  
  roc_final_int <- roc(y_train_all, prob_final_int[, "1"], quiet = TRUE, direction = "<")
  coords_final <- coords(
    roc_final_int, "best",
    ret = c("threshold"),
    transpose = FALSE, best.method = "youden"
  )
  if (is.data.frame(coords_final) && nrow(coords_final) > 1) coords_final <- coords_final[1, ]
  unified_threshold <- coords_final$threshold
  
  cm_int <- confusionMatrix(
    factor(ifelse(preds_cv_int > unified_threshold, 1, 0), levels = c(0, 1)),
    factor(y_train_all, levels = c(0, 1)),
    positive = "1"
  )
  
  int_res <- list(
    roc = roc_int,
    auc = as.numeric(auc(roc_int)),
    ci_lower = as.numeric(roc_int$ci[1]),
    ci_upper = as.numeric(roc_int$ci[3]),
    threshold = unified_threshold,
    sensitivity = as.numeric(cm_int$byClass["Sensitivity"]),
    specificity = as.numeric(cm_int$byClass["Specificity"]),
    accuracy = as.numeric(cm_int$overall["Accuracy"]),
    conf_mat = cm_int,
    n_matched = n_matched,
    k_cv = k_cv
  )
  
  ext_res <- NULL
  
  train_age_min <- min(final_train_matched$Age, na.rm = TRUE)
  train_age_max <- max(final_train_matched$Age, na.rm = TRUE)
  
  sub_test_supported <- apply_train_support_range(sub_test_raw, train_age_min, train_age_max)
  n_test <- nrow(sub_test_supported)
  
  if (n_test > 0) {
    y_test <- sub_test_supported$label
  } else {
    y_test <- c()
  }
  
  if (n_test > 0 && length(unique(y_test)) == 2) {
    final_test_imp <- sub_test_supported
    
    for (col in numeric_cols) {
      if (any(is.na(final_test_imp[[col]]))) final_test_imp[[col]][is.na(final_test_imp[[col]])] <- final_medians[[col]]
    }
    
    pred_ext_only <- tryCatch(
      predict(final_svm_model, newdata = final_test_imp[, final_sel_feat, drop = FALSE], probability = TRUE),
      error = function(e) NULL
    )
    
    if (!is.null(pred_ext_only)) {
      prob_ext <- attr(pred_ext_only, "probabilities")
      if (!is.null(prob_ext) && "1" %in% colnames(prob_ext)) {
        preds_ext_only <- prob_ext[, "1"]
        
        roc_ext <- roc(y_test, preds_ext_only, quiet = TRUE, ci = TRUE, direction = "<")
        ext_threshold <- unified_threshold
        
        cm_ext <- confusionMatrix(
          factor(ifelse(preds_ext_only > ext_threshold, 1, 0), levels = c(0, 1)),
          factor(y_test, levels = c(0, 1)),
          positive = "1"
        )
        
        ext_res <- list(
          roc = roc_ext,
          auc = as.numeric(auc(roc_ext)),
          ci_lower = as.numeric(roc_ext$ci[1]),
          ci_upper = as.numeric(roc_ext$ci[3]),
          threshold = ext_threshold,
          sensitivity = as.numeric(cm_ext$byClass["Sensitivity"]),
          specificity = as.numeric(cm_ext$byClass["Specificity"]),
          accuracy = as.numeric(cm_ext$overall["Accuracy"]),
          conf_mat = cm_ext,
          n_test = n_test,
          train_age_min = train_age_min,
          train_age_max = train_age_max
        )
      }
    }
  }
  
  full_train_imp <- d_train
  full_test_imp  <- d_test
  
  for (col in numeric_cols) {
    if (any(is.na(full_train_imp[[col]]))) full_train_imp[[col]][is.na(full_train_imp[[col]])] <- final_medians[[col]]
    if (nrow(full_test_imp) > 0 && any(is.na(full_test_imp[[col]]))) {
      full_test_imp[[col]][is.na(full_test_imp[[col]])] <- final_medians[[col]]
    }
  }
  
  pred_full_train <- predict(final_svm_model, newdata = full_train_imp[, final_sel_feat, drop = FALSE], probability = TRUE)
  full_preds_train <- attr(pred_full_train, "probabilities")[, "1"]
  
  full_preds_test <- NULL
  if (nrow(full_test_imp) > 0) {
    pred_full_test <- predict(final_svm_model, newdata = full_test_imp[, final_sel_feat, drop = FALSE], probability = TRUE)
    full_preds_test <- attr(pred_full_test, "probabilities")[, "1"]
  }
  
  return(list(
    int = int_res,
    ext = ext_res,
    full_preds_train = full_preds_train,
    full_preds_test = full_preds_test,
    reason = NA,
    seed_records = if(record_seeds) seed_records else NULL
  ))
}

############################################################
## 5. Train disease-specific classifiers
############################################################
diseases <- c("MS", "NMOSD", "MOGAD", "AE")

all_seed_info <- data.frame(
  Disease = character(), Step = character(), Seed_Value = integer(), Description = character(),
  stringsAsFactors = FALSE
)

risk_scores_train_df <- df_train[, c("Age", "Sex", "Diagnosis", "Site_ZZZ")]
risk_scores_test_df  <- df_test[, c("Age", "Sex", "Diagnosis", "Site_ZZZ")]

results_list <- lapply(diseases, function(d) {
  cat("  Analyzing:", d, "vs Rest\n")
  disease_seed <- sample(1:1000000, 1)
  cat("  Random seed:", disease_seed, "\n")
  
  all_seed_info <<- rbind(all_seed_info, data.frame(
    Disease = d, Step = "Global_Seed", Seed_Value = disease_seed,
    Description = paste("Main random seed for", d, "vs Rest"), stringsAsFactors = FALSE
  ))
  
  result <- tryCatch(
    run_disease_classifier(
      df_train, df_test,
      d, use_seed = disease_seed, record_seeds = TRUE
    ),
    error = function(e) {
      cat("  Error occurred but execution continues:", e$message, "\n")
      list(
        int = NULL, ext = NULL, reason = e$message,
        full_preds_train = NULL, full_preds_test = NULL, seed_records = NULL
      )
    }
  )
  
  if (!is.null(result$seed_records) && length(result$seed_records) > 0) {
    for (step_name in names(result$seed_records)) {
      all_seed_info <<- rbind(all_seed_info, data.frame(
        Disease = d, Step = step_name, Seed_Value = as.integer(result$seed_records[[step_name]]),
        Description = paste("Seed for", d, step_name), stringsAsFactors = FALSE
      ))
    }
  }
  
  if (!is.null(result$full_preds_train)) risk_scores_train_df[[paste0(d, "_RiskScore")]] <<- result$full_preds_train
  if (!is.null(result$full_preds_test))  risk_scores_test_df[[paste0(d, "_RiskScore")]]  <<- result$full_preds_test
  
  return(result)
})
names(results_list) <- diseases

res_int_plot <- list()
res_ext_plot <- list()
for (d in names(results_list)) {
  if (!is.null(results_list[[d]]$int) && !is.null(results_list[[d]]$int$roc)) res_int_plot[[d]] <- results_list[[d]]$int
  if (!is.null(results_list[[d]]$ext) && !is.null(results_list[[d]]$ext$roc)) res_ext_plot[[d]] <- results_list[[d]]$ext
}

############################################################
## 6. Define plotting functions
############################################################
plot_roc_curves <- function(res_list, title_text) {
  if (length(res_list) == 0) return()
  colors <- c("#8b60b5", "#344aaa", "#0094d4", "#ad1a71","#e69c24",  "#3CB371", "#FFA07A", "#20B2AA", "#9370DB")
  old_par <- par(lwd = 2, cex.lab = 2, cex.axis = 1.5)
  
  plot(res_list[[1]]$roc, col = colors[1], lwd = 4.5, main = title_text,
       legacy.axes = TRUE, xlab = "1 - Specificity", ylab = "Sensitivity",
       identity = TRUE, identity.lty = 2, identity.lwd = 2)
  
  if (length(res_list) > 1) {
    for (i in 2:length(res_list)) {
      plot(res_list[[i]]$roc, col = colors[i], lwd = 4.5, add = TRUE)
    }
  }
  
  legend_text <- sapply(seq_along(res_list), function(i) {
    d <- names(res_list)[i]
    paste0(
      d, ": AUC=", format(round(res_list[[i]]$auc, 3), nsmall = 3),
      " (", format(round(res_list[[i]]$ci_lower, 3), nsmall = 3),
      "-", format(round(res_list[[i]]$ci_upper, 3), nsmall = 3), ")"
    )
  })
  legend("bottomright", legend = legend_text, col = colors[1:length(res_list)], lwd = 3, cex = 0.9)
  par(old_par)
}



############################################################
## 7. Configure output directory and generate plots
############################################################
out_dir <- "D:/output_folder"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (length(res_int_plot) > 0) {
  pdf(file.path(out_dir, "ROC_Internal_Target_vs_Rest.pdf"), width = 7, height = 7)
  plot_roc_curves(res_int_plot, "ROC Curves Internal (Target vs Rest)")
  dev.off()
  
  png(file.path(out_dir, "ROC_Internal_Target_vs_Rest.png"), width = 7, height = 7, units = "in", res = 300)
  plot_roc_curves(res_int_plot, "ROC Curves Internal (Target vs Rest)")
  dev.off()
}

if (length(res_ext_plot) > 0) {
  pdf(file.path(out_dir, "ROC_External_Target_vs_Rest.pdf"), width = 7, height = 7)
  plot_roc_curves(res_ext_plot, "ROC Curves External (Target vs Rest)")
  dev.off()
  
  png(file.path(out_dir, "ROC_External_Target_vs_Rest.png"), width = 7, height = 7, units = "in", res = 300)
  plot_roc_curves(res_ext_plot, "ROC Curves External (Target vs Rest)")
  dev.off()
}

write.csv(risk_scores_train_df, file.path(out_dir, "Predicted_Risk_Scores_Internal.csv"), row.names = FALSE)
write.csv(risk_scores_test_df, file.path(out_dir, "Predicted_Risk_Scores_External.csv"), row.names = FALSE)

############################################################
## 8. Generate Excel summary with metrics and failure reasons
############################################################
create_excel_summary <- function(res_list, seed_info, output_dir) {
  sum_int <- data.frame(
    Disease = character(),
    AUC = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    Threshold = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    Accuracy = numeric(),
    N_Matched = integer(),
    K_CV = integer(),
    Status = character(),
    Reason = character(),
    stringsAsFactors = FALSE
  )
  
  sum_ext <- data.frame(
    Disease = character(),
    AUC = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    Threshold = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    Accuracy = numeric(),
    N_Test = integer(),
    Train_Age_Min = numeric(),
    Train_Age_Max = numeric(),
    Status = character(),
    Reason = character(),
    stringsAsFactors = FALSE
  )
  
  for (d in names(res_list)) {
    res <- res_list[[d]]
    
    # Internal results
    if (!is.null(res$int)) {
      i_res <- res$int
      sum_int <- rbind(sum_int, data.frame(
        Disease = d,
        AUC = round(i_res$auc, 3),
        CI_Lower = round(i_res$ci_lower, 3),
        CI_Upper = round(i_res$ci_upper, 3),
        Threshold = round(i_res$threshold, 3),
        Sensitivity = round(i_res$sensitivity, 3),
        Specificity = round(i_res$specificity, 3),
        Accuracy = round(i_res$accuracy, 3),
        N_Matched = i_res$n_matched,
        K_CV = i_res$k_cv,
        Status = "Success",
        Reason = "",
        stringsAsFactors = FALSE
      ))
    } else {
      sum_int <- rbind(sum_int, data.frame(
        Disease = d,
        AUC = NA, CI_Lower = NA, CI_Upper = NA,
        Threshold = NA, Sensitivity = NA, Specificity = NA, Accuracy = NA,
        N_Matched = NA, K_CV = NA,
        Status = "Failed",
        Reason = ifelse(is.null(res$reason), "Unknown failure", res$reason),
        stringsAsFactors = FALSE
      ))
    }
    
    # External results
    if (!is.null(res$ext)) {
      e_res <- res$ext
      sum_ext <- rbind(sum_ext, data.frame(
        Disease = d,
        AUC = round(e_res$auc, 3),
        CI_Lower = round(e_res$ci_lower, 3),
        CI_Upper = round(e_res$ci_upper, 3),
        Threshold = round(e_res$threshold, 3),
        Sensitivity = round(e_res$sensitivity, 3),
        Specificity = round(e_res$specificity, 3),
        Accuracy = round(e_res$accuracy, 3),
        N_Test = e_res$n_test,
        Train_Age_Min = round(e_res$train_age_min, 3),
        Train_Age_Max = round(e_res$train_age_max, 3),
        Status = "Success",
        Reason = "",
        stringsAsFactors = FALSE
      ))
    } else {
      sum_ext <- rbind(sum_ext, data.frame(
        Disease = d,
        AUC = NA, CI_Lower = NA, CI_Upper = NA,
        Threshold = NA, Sensitivity = NA, Specificity = NA, Accuracy = NA,
        N_Test = NA, Train_Age_Min = NA, Train_Age_Max = NA,
        Status = "Failed/Not available",
        Reason = ifelse(is.null(res$reason), "External result unavailable", res$reason),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  sum_int <- sum_int[order(sum_int$AUC, decreasing = TRUE, na.last = TRUE), ]
  sum_ext <- sum_ext[order(sum_ext$AUC, decreasing = TRUE, na.last = TRUE), ]
  
  wb <- createWorkbook()
  addWorksheet(wb, "Metrics_Internal")
  writeData(wb, "Metrics_Internal", sum_int)
  
  addWorksheet(wb, "Metrics_External")
  writeData(wb, "Metrics_External", sum_ext)
  
  addWorksheet(wb, "Random_Seeds")
  writeData(wb, "Random_Seeds", seed_info)
  
  saveWorkbook(wb, file.path(output_dir, "ROC_results_summary_Int_and_Ext.xlsx"), overwrite = TRUE)
}

create_excel_summary(results_list, all_seed_info, out_dir)


############################################################
## 9. Print analysis summary
############################################################
cat("Independent internal and external validation completed (Target vs Rest).\n")
cat("External test site(s):", paste(test_sites, collapse = ", "), "\n")
cat("Output directory:", out_dir, "\n")
cat("Generated key files:\n")
cat("  1. ROC_Internal / External_Target_vs_Rest: ROC curves for successful target-vs-rest classifiers\n")
cat("  2. Predicted_Risk_Scores_Internal / External.csv: risk score tables for internal and external sets\n")
cat("  3. ROC_results_summary_Int_and_Ext.xlsx: metrics, success/failure status, and failure reasons\n")