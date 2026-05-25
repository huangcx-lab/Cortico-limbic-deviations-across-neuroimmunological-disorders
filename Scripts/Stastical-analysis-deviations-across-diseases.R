# =====================================================================
# 1. Environment setup and package loading
# =====================================================================
rm(list=ls())


library(ggplot2)
library(lmerTest)
library(effectsize)
library(ggsci)
library(ggdist)
library(MatchIt)
library(dplyr)


input_file <- "D:/input_folder/data.csv"
output_dir <- "D:/output_folder"


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

# =====================================================================
# 2. Data loading and preparation
# =====================================================================
cat("Loading data, please wait...\n")
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# combine LGI1, NMDA, GAD65, GABA, and CASPR2 into the AE group
ae_subtypes <- c('LGI1', 'NMDA', 'GAD65', 'GABA', 'CASPR2')
ae_data <- data[data$Diagnosis %in% ae_subtypes, ]
if (nrow(ae_data) > 0) {
  ae_data$Diagnosis <- 'AE'
  data <- rbind(data, ae_data)  
}
# ------------------------------------------------------------------

data$Sex <- as.factor(data$Sex)
data$Site_ZZZ <- as.factor(data$Site_ZZZ)
data$Diagnosis <- as.factor(data$Diagnosis)


disease_groups <- c('MS', 'AQP4Pos_NMOSD', 'MOGAD', 'AE', 'LGI1', 'NMDA')

features <- grep("^Quant_", colnames(data), value = TRUE)
cat(sprintf("Successfully identified %d feature variables.\n", length(features)))

all_results_list <- list()
row_idx <- 1

# =====================================================================
# 3. Core analysis loop
# =====================================================================
for (i_feature in features) {
  cat(sprintf("\n>>> Processing feature: %s", i_feature))
  
  for (disease in disease_groups) {
    data_disease <- data[data$Diagnosis == disease & !is.na(data$Diagnosis), ]
    if(nrow(data_disease) < 3) next 
    
    min_age <- min(data_disease$Age, na.rm = TRUE)
    max_age <- max(data_disease$Age, na.rm = TRUE)
    tem_data1 <- data[data$Diagnosis %in% c('HC', disease), ]
    tem_data1 <- tem_data1[tem_data1$Age >= min_age & tem_data1$Age <= max_age, ]
    tem_data1 <- tem_data1[!is.na(tem_data1[[i_feature]]) & !is.na(tem_data1$Age) & !is.na(tem_data1$Sex), ]
    
    if(nrow(tem_data1[tem_data1$Diagnosis == 'HC',]) < 5) next
    
    tem_data1$Diagnosis <- factor(tem_data1$Diagnosis, levels = c('HC', disease))
    match_success <- tryCatch({
      match.it <- matchit(Diagnosis ~ Age + Sex + Site_ZZZ, data = tem_data1, method = "nearest", ratio = 1)
      tem_matched <- match.data(match.it)
      TRUE
    }, error = function(e) FALSE)
    
    if(!match_success) next
    
    G1 <- tem_matched[[i_feature]][tem_matched$Diagnosis == 'HC']
    G2 <- tem_matched[[i_feature]][tem_matched$Diagnosis == disease]
    
    if(length(G1) < 3 || length(G2) < 3 || var(G1, na.rm=T) == 0 || var(G2, na.rm=T) == 0) next
    
    N_HC <- length(G1)
    N_Disease <- length(G2)
    
    form <- as.formula(paste0(i_feature, " ~ Diagnosis"))
    tem_stat <- lm(form, data = tem_matched)
    tem_stats <- summary(tem_stat)
    
    tem_d <- cohens_d(form, data = tem_matched)
    tem_g <- hedges_g(form, data = tem_matched)
    tem_delta <- glass_delta(form, data = tem_matched)
    tem_eta2 <- eta_squared(aov(form, data = tem_matched))
    
    tem_matched$label <- ifelse(tem_matched$Diagnosis == disease, 1, 0)
    log_form <- as.formula(paste0("label ~ ", i_feature))
    tem_logs <- tryCatch(glm(log_form, family = binomial, data = tem_matched), error = function(e) NULL)
    
    or_val <- NA; or_low <- NA; or_high <- NA
    if(!is.null(tem_logs)){
      or_val <- tem_logs$coefficients[2]
      ci_log <- tryCatch(suppressMessages(confint(tem_logs)), error = function(e) matrix(NA, 2, 2))
      or_low <- ci_log[2,1]; or_high <- ci_log[2,2]
    }
    
    tem_op <- tryCatch(p_overlap(G1, G2, parametric = FALSE), 
                       error = function(e) list(Overlap=NA, CI_low=NA, CI_high=NA))
    
    tem_rb <- rank_biserial(G1, G2)
    tem_cld <- cliffs_delta(G1, G2)
    tem_u1 <- cohens_u1(G1, G2)
    tem_wilcox <- wilcox.test(G1, G2)
    
    res_row <- data.frame(
      feature = i_feature, 
      disease = disease,
      N_HC = N_HC, 
      N_Disease = N_Disease, 
      pvalue_raw = tem_stats$coefficients[2,4], 
      tvalue = tem_stats$coefficients[2,3],
      degree = tem_stats$df[2], 
      Fvalue = ifelse(is.null(tem_stats$fstatistic), NA, tem_stats$fstatistic[1]),
      Adjust_R2 = tem_stats$adj.r.squared,
      cohens_d = tem_d$Cohens_d, cohens_d_CI_low = tem_d$CI_low, cohens_d_CI_high = tem_d$CI_high,
      hedges_g = tem_g$Hedges_g, hedges_g_CI_low = tem_g$CI_low, hedges_g_CI_high = tem_g$CI_high,
      glass_delta = tem_delta$Glass_delta, glass_delta_CI_low = tem_delta$CI_low, glass_delta_CI_high = tem_delta$CI_high,
      Eta2 = tem_eta2$Eta2, Eta2_CI_low = tem_eta2$CI_low, Eta2_CI_high = tem_eta2$CI_high,
      OR_LogOdds = or_val, OR_CI_low = or_low, OR_CI_high = or_high,
      Overlap = tem_op$Overlap, Overlap_CI_low = tem_op$CI_low, Overlap_CI_high = tem_op$CI_high,
      rank_biserial = tem_rb$r_rank_biserial, rank_biserial_CI_low = tem_rb$CI_low, rank_biserial_CI_high = tem_rb$CI_high,
      cliffs_delta = tem_cld$r_rank_biserial, cliffs_delta_CI_low = tem_cld$CI_low, cliffs_delta_CI_high = tem_cld$CI_high,
      cohens_u1 = tem_u1$Cohens_U1, cohens_u1_CI_low = tem_u1$CI_low, cohens_u1_CI_high = tem_u1$CI_high,
      wilcox_test_W = tem_wilcox$statistic, 
      wilcox_p_raw = tem_wilcox$p.value
    )
    
    all_results_list[[row_idx]] <- res_row
    row_idx <- row_idx + 1
  }
  

}

# =====================================================================
# 4. Multiple-comparison correction (FDR) and result formatting
# =====================================================================
final_df <- do.call(rbind, all_results_list)

# Define significance-star labeling function
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05)  return("*")
  return("")
}

# Apply FDR correction separately within each disease group across brain regions
final_df <- final_df %>%
  group_by(disease) %>%
  mutate(
    # Parametric-test correction
    p_FDR_LM = p.adjust(pvalue_raw, method = "fdr"),
    sig_FDR_LM = sapply(p_FDR_LM, get_stars),
    
    # Non-parametric-test correction
    p_FDR_Wilcox = p.adjust(wilcox_p_raw, method = "fdr"),
    sig_FDR_Wilcox = sapply(p_FDR_Wilcox, get_stars)
  ) %>%
  ungroup()


final_df <- final_df %>%
  relocate(p_FDR_LM, sig_FDR_LM, .after = pvalue_raw) %>%
  relocate(p_FDR_Wilcox, sig_FDR_Wilcox, .after = wilcox_p_raw)

# Save final results table
write.csv(final_df, "All_Features_Stats_Results_FDR.csv", row.names = FALSE)

cat("\n\nAnalysis complete! Statistical results have been saved to: All_Features_Stats_Results_FDR.csv\n")
cat("FDR correction has been performed within each disease group, with significance stars added.\n")