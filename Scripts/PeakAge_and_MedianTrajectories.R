# -----------------------------
# 1. Install and load packages
# -----------------------------
required_packages <- c(
  "ggplot2", "dplyr", "stringr", "readr", "forcats",
  "tidyr", "purrr", "tools", "scales", "ggpattern"
)

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(forcats)
library(tidyr)
library(purrr)
library(tools)
library(scales)
library(ggpattern)

# -----------------------------
# 2. Set input and output paths
# -----------------------------
input_dir <- "D:/input_folder"
output_dir <- "D:/output_folder"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

median_raw_dir <- file.path(output_dir, "Median_Raw_Curves")
if (!dir.exists(median_raw_dir)) {
  dir.create(median_raw_dir, recursive = TRUE)
}

# -----------------------------
# 3. Read all CSV files
# -----------------------------
csv_files <- list.files(
  path = input_dir,
  pattern = "_p2\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("No matching CSV files were found in the input directory. Please check the input path.")
}

# -----------------------------
# 4. Define brain region order
# -----------------------------
target_order <- c(
  "Hippocampus",
  "Amygdala",
  "Accumbens",
  "anteriorcingulate",
  "posteriorcingulate",
  "entorhinal",
  "orbitofrontal"
)


normalize_region_name <- function(region_base_raw) {
  x <- region_base_raw
  x_low <- tolower(x)
  
  if (str_detect(x_low, "hippo")) return("Hippocampus")
  if (str_detect(x_low, "amyg")) return("Amygdala")
  if (str_detect(x_low, "accumb")) return("Accumbens")
  if (str_detect(x_low, "anteriorcing")) return("anteriorcingulate")
  if (str_detect(x_low, "posteriorcing")) return("posteriorcingulate")
  if (str_detect(x_low, "entorh")) return("entorhinal")
  if (str_detect(x_low, "orbitofr")) return("orbitofrontal")
  
  return(region_base_raw)
}

# -----------------------------
# 5. Parse file information
# -----------------------------
parse_file_info <- function(file_path) {
  file_name <- basename(file_path)
  file_name_no_ext <- file_path_sans_ext(file_name)
  
  region_full_raw <- str_remove(file_name_no_ext, "_数据框_p2$")
  
  hemisphere <- case_when(
    str_detect(region_full_raw, "^lh[_\\.-]") ~ "lh",
    str_detect(region_full_raw, "^rh[_\\.-]") ~ "rh",
    TRUE ~ "unknown"
  )
  
  region_base_raw <- region_full_raw %>%
    str_remove("^lh[_\\.-]") %>%
    str_remove("^rh[_\\.-]")
  
  region_group <- normalize_region_name(region_base_raw)
  
  hemisphere_cn <- case_when(
    hemisphere == "lh" ~ "LH",
    hemisphere == "rh" ~ "RH",
    TRUE ~ "UNK"
  )
  
  display_label <- paste0(region_group, "_", hemisphere_cn)
  
  list(
    file_name = file_name,
    region_full_raw = region_full_raw,
    hemisphere = hemisphere,
    region_base_raw = region_base_raw,
    region_group = region_group,
    display_label = display_label
  )
}

# -----------------------------
# 6. Define colors by brain region
# -----------------------------
region_color_map <- c(
  "Hippocampus"        = "#FFD340",
  "Amygdala"           = "#c25109",
  "Accumbens"          = "#807381",
  "anteriorcingulate"  = "#0094e1",
  "posteriorcingulate" = "#394cbb",
  "entorhinal"         = "#bf0a77",
  "orbitofrontal"      = "#750A87"
)


all_region_groups_detected <- unique(map_chr(csv_files, ~ parse_file_info(.x)$region_group))
unknown_regions <- setdiff(all_region_groups_detected, names(region_color_map))

if (length(unknown_regions) > 0) {
  extra_cols <- hue_pal()(length(unknown_regions))
  names(extra_cols) <- unknown_regions
  region_color_map <- c(region_color_map, extra_cols)
}

# -----------------------------
# 7. Read data and extract peak age
# -----------------------------
summary_list <- list()
curve_list <- list()

for (i in seq_along(csv_files)) {
  file <- csv_files[i]
  info <- parse_file_info(file)
  
  df <- tryCatch(
    read_csv(file, show_col_types = FALSE),
    error = function(e) {
      message("Failed to read file:", basename(file))
      return(NULL)
    }
  )
  
  if (is.null(df)) next
  
  if (!all(c("Age", "median") %in% colnames(df))) {
    message("Skipped file because required columns Age or median are missing:", basename(file))
    next
  }
  
  df2 <- df %>%
    select(Age, median) %>%
    filter(!is.na(Age), !is.na(median)) %>%
    arrange(Age)
  
  if (nrow(df2) == 0) {
    message("Skipped file because no valid Age/median data were available:", basename(file))
    next
  }
  

  peak_row <- df2 %>%
    filter(median == max(median, na.rm = TRUE)) %>%
    arrange(Age) %>%
    slice(1)
  
  peakage <- peak_row$Age[1]
  maxmedian <- peak_row$median[1]
  
  summary_list[[length(summary_list) + 1]] <- data.frame(
    file_name = info$file_name,
    region_full_raw = info$region_full_raw,
    region_base_raw = info$region_base_raw,
    region_group = info$region_group,
    hemisphere = info$hemisphere,
    display_label = info$display_label,
    peakage = peakage,
    maxmedian = maxmedian,
    stringsAsFactors = FALSE
  )
  
  curve_list[[length(curve_list) + 1]] <- df2 %>%
    mutate(
      file_name = info$file_name,
      region_group = info$region_group,
      hemisphere = info$hemisphere,
      display_label = info$display_label
    )
}

summary_df <- bind_rows(summary_list)
curve_df <- bind_rows(curve_list)

if (nrow(summary_df) == 0) {
  stop("No files were successfully processed. Please check the file format.")
}

# -----------------------------
# 8. Sort data
# -----------------------------
summary_df <- summary_df %>%
  mutate(
    region_group = factor(region_group, levels = target_order),
    hemi_order = case_when(
      hemisphere == "lh" ~ 1,
      hemisphere == "rh" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(region_group, hemi_order)

curve_df <- curve_df %>%
  mutate(
    region_group = factor(region_group, levels = target_order),
    hemi_order = case_when(
      hemisphere == "lh" ~ 1,
      hemisphere == "rh" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(region_group, hemi_order, Age)


curve_df <- curve_df %>%
  group_by(display_label) %>%
  mutate(
    median_min = min(median, na.rm = TRUE),
    median_max = max(median, na.rm = TRUE),
    percentage = ifelse(
      median_max > median_min,
      (median - median_min) / (median_max - median_min) * 100,
      0
    ),
    linetype_hemi = ifelse(hemisphere == "lh", "solid", "dashed")
  ) %>%
  ungroup() %>%
  select(-median_min, -median_max)

# -----------------------------
# 9. Create color mappings
# -----------------------------
summary_df <- summary_df %>%
  mutate(
    fill_color = region_color_map[as.character(region_group)]
  )

curve_df <- curve_df %>%
  mutate(
    line_color = region_color_map[as.character(region_group)]
  )


fill_color_map <- summary_df$fill_color
names(fill_color_map) <- summary_df$display_label

line_color_map <- curve_df %>%
  distinct(display_label, line_color) %>%
  { setNames(.$line_color, .$display_label) }

# -----------------------------
# 10. Export summary table
# -----------------------------
summary_outfile <- file.path(output_dir, "p2_peakage_summary_enhanced.csv")
write_csv(summary_df %>% arrange(region_group, hemi_order), summary_outfile)


# -----------------------------
# 11. Create horizontal peak-age bar plot
# -----------------------------

bar_width <- 1
inner_gap <- 0.1
between_gap <- 0.7

summary_df_bar <- summary_df %>%
  mutate(
    region_group = factor(region_group, levels = target_order),
    hemi_order_bar = ifelse(hemisphere == "lh", 1, 2)
  ) %>%
  arrange(region_group, hemi_order_bar) %>%
  group_by(region_group) %>%
  mutate(within_group_order = row_number()) %>%
  ungroup() %>%
  arrange(region_group, within_group_order)

pos_vec <- numeric(nrow(summary_df_bar))
current_pos <- 1
for (i in 1:nrow(summary_df_bar)) {
  pos_vec[i] <- current_pos
  if (i < nrow(summary_df_bar)) {
    if (summary_df_bar$region_group[i] == summary_df_bar$region_group[i + 1]) {
      current_pos <- current_pos + bar_width + inner_gap
    } else {
      current_pos <- current_pos + bar_width + between_gap
    }
  }
}
summary_df_bar$pos_bar <- pos_vec

summary_df_bar <- summary_df_bar %>%
  mutate(pos_bar_rev = max(pos_bar) - pos_bar + min(pos_bar))

summary_df_bar <- summary_df_bar %>%
  mutate(
    label_text = paste0(round(peakage, 2), "yr / ", round(maxmedian, 2), "mm³")
  )

summary_df_bar <- summary_df_bar %>%
  mutate(
    region_name_clean = case_when(
      region_group == "anteriorcingulate" ~ "Anteriorcingulate",
      region_group == "posteriorcingulate" ~ "Posteriorcingulate",
      region_group == "entorhinal" ~ "Entorhinal",
      region_group == "orbitofrontal" ~ "Orbitofrontal",
      TRUE ~ as.character(region_group)
    ),
    hemi_clean = ifelse(hemisphere == "lh", "L", "R"),
    display_label_clean = paste0(region_name_clean, "_", hemi_clean)
  )

summary_df_bar <- summary_df_bar %>%
  mutate(
    border_color = ifelse(hemisphere == "lh", fill_color, "black")
  )

bar_breaks_rev <- summary_df_bar$pos_bar_rev
bar_labels <- summary_df_bar$display_label_clean

xmax_peak <- max(summary_df_bar$peakage, na.rm = TRUE)
label_offset <- xmax_peak * 0.03
x_limit_peak <- xmax_peak * 1.5

p_bar_peak <- ggplot(summary_df_bar, aes(x = pos_bar_rev, y = peakage)) +
  ggpattern::geom_col_pattern(
    aes(fill = region_group,
        pattern = hemisphere,
        color = border_color),
    pattern_fill = "white",
    pattern_color = NA,
    pattern_density = 0.3,
    pattern_spacing = 0.02,
    pattern_angle = 45,
    pattern_linetype = "dashed",
    width = bar_width,
    linewidth = 0.5
  ) +
  scale_fill_manual(values = region_color_map) +
  scale_pattern_manual(values = c("lh" = "stripe", "rh" = "none")) +
  scale_color_identity() +
  geom_text(
    aes(y = peakage + label_offset, label = label_text),
    hjust = 0,
    size = 6
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, x_limit_peak), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = bar_breaks_rev,
    labels = bar_labels,
    expand = expansion(add = c(0.3, 0.3))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Peak Age (yr)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    legend.position = "right",
    guides(pattern = "none")
  )

ggsave(
  filename = file.path(output_dir, "01_peakage_horizontal_bar_grouped_LH_RH.png"),
  plot = p_bar_peak,
  width = 12,
  height = max(6, 0.85 * nrow(summary_df_bar)),
  dpi = 300
)



# -----------------------------
# 12. Plot overlaid normalized median trajectories
# -----------------------------
curve_df_overlay <- curve_df %>%
  filter(!is.na(Age), Age > 0)

p_percent_overlay <- ggplot(curve_df_overlay, aes(x = Age, y = percentage, 
                                                  color = display_label, linetype = linetype_hemi)) +

  geom_line(linewidth = 2, alpha = 0.95, lineend = "round") +
  scale_color_manual(values = fill_color_map) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed"),
                        labels = c("LH", "RH"), name = "Hemisphere") +
  scale_x_continuous(
    trans = "sqrt",
    breaks = c(6, 18, 35, 80),
    labels = c("6", "18", "35", "80"),
    minor_breaks = NULL
  ) +
  scale_y_continuous(limits = c(-5, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(
    title = "All Age-Percentage Curves Overlay (P2)",
    subtitle = "Percentage = (median - min) / (max - min) * 100",
    x = "Age (Years)",
    y = "Normalized Trajectory (% of range)",
    color = "Region"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.text = element_text(size = 14)
  )

ggsave(
  filename = file.path(output_dir, "10_percentage_curves_overlay.png"),
  plot = p_percent_overlay,
  width = 12,
  height = 8,
  dpi = 300
)


# -----------------------------
# 13. Plot faceted normalized median trajectories
# -----------------------------
p_percent_facet <- ggplot(curve_df, aes(x = Age, y = percentage, 
                                        color = display_label, linetype = linetype_hemi)) +
  geom_line(linewidth = 1.2, alpha = 0.95, lineend = "round") +
  scale_color_manual(values = fill_color_map) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed"),
                        labels = c("LH", "RH"), name = "Hemisphere") +

  scale_x_continuous(
    trans = "sqrt",
    breaks = c(6, 18, 35, 80),
    labels = c("6", "18", "35", "80"),
    minor_breaks = NULL
  ) +
  scale_y_continuous(limits = c(-5, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  facet_wrap(~ region_group, scales = "fixed", ncol = 3) +
  labs(
    title = "Age-Percentage Curves by Brain Region (P2)",
    subtitle = "Percentage = (median - min) / (max - min) * 100; LH solid, RH dashed",
    x = "Age (Years)",
    y = "Normalized Trajectory (% of range)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(output_dir, "11_percentage_curves_facet_by_region.png"),
  plot = p_percent_facet,
  width = 13,
  height = 9,
  dpi = 300
)


# -----------------------------
# 14. Plot median curves by brain region
# -----------------------------
for (rg in target_order) {
  d_sub <- curve_df %>%
    filter(region_group == rg)
  
  if (nrow(d_sub) == 0) next
  
  p_median_single <- ggplot(d_sub, aes(x = Age, y = median, 
                                       color = display_label, linetype = linetype_hemi)) +
    geom_line(linewidth = 1.5, alpha = 0.95, lineend = "round") +
    scale_color_manual(values = fill_color_map) +
    scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed"),
                          labels = c("LH", "RH"), name = "Hemisphere") +
    scale_x_continuous(
      trans = "sqrt",
      breaks = c(6, 18, 35, 80),
      labels = c("6", "18", "35", "80")
    ) +
    labs(
      title = paste0("Median raw trajectories - ", rg),
      x = "Age (Years)",
      y = "Median Volume",
      color = "Region"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.2),
      axis.text = element_text(size = 14)
    )
  
  ggsave(
    filename = file.path(median_raw_dir, paste0("20_median_raw_", rg, ".png")),
    plot = p_median_single,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# -----------------------------
# 15. Create faceted median curves
# -----------------------------
p_median_facet <- ggplot(curve_df, aes(x = Age, y = median, 
                                       color = display_label, linetype = linetype_hemi)) +
  geom_line(linewidth = 1.2, alpha = 0.95, lineend = "round") +
  scale_color_manual(values = fill_color_map) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "33"),
                        labels = c("LH", "RH"), name = "Hemisphere") +
  scale_x_continuous(
    trans = "sqrt",
    breaks = c(6, 18, 35, 80),
    labels = c("6", "18", "35", "80")
  ) +
  facet_wrap(~ region_group, scales = "free", ncol = 3, labeller = as_labeller(stringr::str_to_title)) +
  labs(
    title = "Median raw trajectories by brain region",
    subtitle = "LH solid, RH dashed",
    x = "Age (Years)",
    y = "Median Volume"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(t = 6, b = 6)),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(median_raw_dir, "21_median_raw_facet_by_region.png"),
  plot = p_median_facet,
  width = 13,
  height = 9,
  dpi = 300
)


report_file <- file.path(output_dir, "README_plot_outputs.txt")
report_text <- c("Finish") 
writeLines(report_text, con = report_file)


cat("Finish!\n")