#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
#### PER SITE STATISTICS FOR INDIAN JARVIS LOW-COVERAGE SEQUENCING
################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#-------------------------------------------------------------------------------
#### LOAD LIBRARIES
#-------------------------------------------------------------------------------

library(tidyverse)
library(scales)
library(patchwork)
library(ggtext)     
library(tidyverse)


#-------------------------------------------------------------------------------
#### SET QC DIRECTORIES
#-------------------------------------------------------------------------------

qc_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/qc_stats/"


#-------------------------------------------------------------------------------
#### FUNCTION FOR UNIFORM THEME
#-------------------------------------------------------------------------------

theme_qc <- function() {
  theme_minimal(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = 11, color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.caption = element_text(color = "grey50", size = 9),
      panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5),  # Add border
      plot.margin = margin(15, 15, 15, 15)
    )
}


#-------------------------------------------------------------------------------
#### Fig 1. DEPTH DISTRIBUTION (FILTERED FOR LOW-COVERAGE REALITY)
#-------------------------------------------------------------------------------

# Read the table
dp <- read_table(paste0(qc_dir, "all_dp_clean.txt"), col_names = "depth") %>%
  filter(depth <= 10)  ## 🔥 NEW: Cap at 10x for 0.8x data

# Plot
p_dp <- ggplot(dp, aes(x = depth)) +
  geom_histogram(binwidth = 1, fill = "#1f77b4", alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  ## 🔥 NEW: Critical threshold
  scale_y_continuous(labels = comma, trans = "log10") +
  labs(
    title = "Variant-level Depth Distribution (0.8x Target)",
    subtitle = "Red line = 1× coverage (minimum for calling)",
    x = "Depth (×)",
    y = "Count (log10)"
  ) +
  theme_qc()


#-------------------------------------------------------------------------------
#### Fig 2. ALLELE FREQUENCY SPECTRUM (FILTERED FOR LOW-AF NOISE)
#-------------------------------------------------------------------------------

# Read the table
af_df <- read_table(paste0(qc_dir, "af_bins_final.txt"), 
                    col_names = c("bin_max", "site_count")) %>%
  filter(bin_max >= 0.01)  ## 🔥 NEW: Remove ultra-rare bins (likely noise)

# Plot
p_af <- ggplot(af_df, aes(x = bin_max, y = site_count)) +
  geom_col(fill = "#ff7f0e", width = 0.7) +
  scale_x_continuous(
    breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1),
    labels = percent_format(accuracy = 1),
    trans = "log10"
  ) +
  labs(
    title = "Minor Allele Frequency Spectrum",
    subtitle = "Excluding AF < 1% (likely sequencing artifacts)",
    x = "Allele Frequency (log scale)",
    y = "Number of Variants"
  ) +
  theme_qc()


#-------------------------------------------------------------------------------
#### Fig 3. MISSINGNESS PER SAMPLE (HIGHLIGHT OUTLIERS)
#-------------------------------------------------------------------------------

# Read the table
miss_df <- read_table(paste0(qc_dir, "missing_counts.tsv"),
                      col_names = c("sample", "missing_gt")) %>%
  mutate(outlier = missing_gt > quantile(missing_gt, 0.95))  ## 🔥 NEW: Flag outliers

# Plot
p_miss <- ggplot(miss_df, aes(reorder(sample, missing_gt), missing_gt, fill = outlier)) +
  geom_col() +
  scale_fill_manual(values = c("#2ca02c", "red"), guide = "none") +
  coord_flip() +
  labs(
    title = "Missing Genotypes per Sample",
    subtitle = "<span style='color:red'>Red</span> = Top 5% outliers",
    x = NULL,
    y = "Missing Genotypes (./.)"
  ) +
  theme_qc() +
  theme(plot.subtitle = element_markdown())  ## 🔥 NEW: Allows HTML in subtitle



#-------------------------------------------------------------------------------
#### Fig 4. Ti/Tv RATIO BY CHROMOSOME (WITH EXPECTED RANGE)
#-------------------------------------------------------------------------------

# Read the table
titv_chr <- read_table(paste0(qc_dir, "titv_per_chr.tsv"),
                       col_names = c("chr", "titv")) %>%
  mutate(chr = factor(chr, levels = 1:10))

# Plot
p_titv <- ggplot(titv_chr, aes(chr, titv, fill = titv)) +
  geom_col() +
  geom_hline(yintercept = c(2.0, 2.2), linetype = "dashed", color = "grey30") +  ## 🔥 NEW: Expected range
  scale_fill_gradient(low = "#d62728", high = "#9467bd") +
  labs(
    title = "Ti/Tv Ratio by Chromosome",
    subtitle = "Dashed lines = Expected range (2.0-2.2 for humans)",
    x = "Chromosome",
    y = "Transition/Transversion Ratio"
  ) +
  theme_qc() +
  theme(legend.position = "none")


#-------------------------------------------------------------------------------
#### COMBINE ALL PLOTS INTO PUBLICATION-READY PANEL
#-------------------------------------------------------------------------------

# Plot
quartz()
(p_dp + p_af) / (p_miss + p_titv) +  ## 🔥 NEW: Uses patchwork
  plot_annotation(
    title = "Low-Coverage (0.8×) Sequencing QC Metrics",
    caption = paste("Data:", basename(qc_dir), "| Filter: DP ≤ 10×, AF ≥ 1%"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save
ggsave("qc_metrics_lowcov.png", width = 14, height = 10, dpi = 300)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
#### PER SAMPLE STATISTICS FOR INDIAN JARVIS LOW-COVERAGE SEQUENCING
################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#-------------------------------------------------------------------------------
#### Load Per-Sample Depth Data
#-------------------------------------------------------------------------------

# Load the table
depth_df <- read_tsv(
  file.path(qc_dir, "per_sample_mean_depth.tsv"),
  col_names = c("sample", "mean_depth", "n_genotypes")
) %>%
  mutate(
    coverage_group = cut(
      mean_depth,
      breaks = c(0, 0.5, 1, 1.5, Inf),
      labels = c("<0.5×", "0.5-1×", "1-1.5×", ">1.5×")
    )
  )

# Calculate global average depth (for reference line)
global_avg_depth <- mean(depth_df$mean_depth)

# Count samples per coverage group (for bar labels)
group_counts <- depth_df %>%
  count(coverage_group) %>%
  mutate(label = paste0(coverage_group, "\n(n=", n, ")"))

group_counts <- depth_df %>%
  group_by(coverage_group) %>%
  summarise(
    n            = n(),
    median_calls = median(n_genotypes),
    .groups      = "drop"
  )

#-------------------------------------------------------------------------------
#### Fig 1.  Depth vs. Genotype Count
#-------------------------------------------------------------------------------

# Calculate ACTUAL Median Genotypes per Coverage Group
depth_stats <- depth_df %>%
  group_by(coverage_group) %>%
  summarise(
    n_samples = n(),
    median_genotypes = median(n_genotypes),
    mean_depth = round(mean(mean_depth), 2)
  ) %>%
  mutate(
    # Format labels with proper scaling
    scatter_label = sprintf(
      "%s (n=%d, median=%s)", 
      coverage_group,
      n_samples,
      scales::comma(median_genotypes, scale = 1e-6, suffix = "M")  # Convert to millions
    ),
    bar_label = sprintf(
      "%s (n=%d, avg=%.2f×)", 
      coverage_group,
      n_samples,
      mean_depth
    )
  )


#-------------------------------------------------------------------------------
#### Fig 2.  Depth vs. Genotyping Attempts
#-------------------------------------------------------------------------------

# Plot
p_scatter <- ggplot(depth_df, aes(mean_depth, n_genotypes / 1e6, color = coverage_group)) +  # Scale Y-axis to millions
  geom_point(alpha = 0.8, size = 3, shape = 16) +
  scale_color_manual(
    values = c("<0.5×" = "dodgerblue", "0.5-1×" = "goldenrod", 
               "1-1.5×" = "purple", ">1.5×" = "olivedrab3"),
    name = "Coverage Group",
    labels = depth_stats$scatter_label  # Use pre-formatted labels
  ) +
  scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
  scale_y_continuous(
    labels = scales::comma,  # Already scaled to millions
    name = "Number of Genotypes Called (Millions)"
  ) +
  labs(
    title = "Depth vs. Genotyping Attempts",
    x = "Mean Depth (×)"
  ) +
  theme_qc() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )


#-------------------------------------------------------------------------------
#### Fig 2.  Ranked Samples
#-------------------------------------------------------------------------------

# Plot
p_ranked <- ggplot(depth_df, aes(reorder(sample, -mean_depth), mean_depth, fill = coverage_group)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.1) +
  # Add group summaries at the top
  geom_text(
    data = depth_df %>% 
      group_by(coverage_group) %>% 
      summarise(
        y_pos = max(mean_depth) * 1.1,
        label = paste0(
          "n=", n(), 
          "\nAvg=", round(mean(mean_depth), 2), "×"
        )
      ),
    aes(x = Inf, y = y_pos, label = label, color = coverage_group),
    hjust = 1.1, size = 3.2, lineheight = 0.8
  ) +
  scale_fill_manual(
    values = c("dodgerblue", "goldenrod", "purple", "olivedrab3"),
    name = "Coverage Group"
  ) +
  scale_color_manual(values = c("dodgerblue4", "goldenrod4", "purple4", "olivedrab4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # More space for text
  labs(
    title = "Sample Depth Ranking",
    x = "Samples (Sorted by Depth)",
    y = "Mean Depth (×)"
  ) +
  theme_qc() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = "none") # Hide color legend for text



#-------------------------------------------------------------------------------
#### COMBINE ALL PLOTS INTO PUBLICATION-READY PANEL
#-------------------------------------------------------------------------------

# Final plot
final_plot <- (p_density + p_scatter) / p_ranked +
  plot_annotation(
    title = "Per-Sample Sequencing Depth Quality Control",
    subtitle = "Low-coverage whole-genome sequencing (0.8× target)",
    caption = paste(
      "Data:", nrow(depth_df), "samples |",
      "Global mean depth:", round(global_avg_depth, 2), "× |",
      "Dashed red line = Target coverage"
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
      plot.caption = element_text(size = 9, hjust = 0.5)
    )
  )

quartz()
final_plot


# Save 
ggsave(
  "Figure3_per_sample_depth_qc.png",
  plot = final_plot,
  width = 10,       # Nature single-column width (inches)
  height = 10,
  dpi = 600,
  bg = "white"
)

getwd()


#-------------------------------------------------------------------------------
#### TABLE FOR STATISTICS
#-------------------------------------------------------------------------------

# Table of statistics
depth_df %>% 
  group_by(coverage_group) %>% 
  summarise(
    n = n(),
    perc = n/nrow(depth_df)*100,
    mean_dp = mean(mean_depth)
  ) %>% 
  mutate(
    label = sprintf("%s: %d samples (%.1f%%) | Avg depth=%.2f×", 
                    coverage_group, n, perc, mean_dp)
  ) %>% 
  pull(label) %>% 
  paste(collapse = "\n") %>% 
  cat()

