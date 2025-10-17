# ──────────────────────────────────────────────
# 0.  Libraries
# ──────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(viridis)     # colour scales
library(patchwork)   # combine the two panels
library(GenWin)      # splineAnalyze()
library(vroom)
library(dplyr)
library(ggplot2)




# ──────────────────────────────────────────────
# 1.  Read & pre‑filter raw table   -------------
# ──────────────────────────────────────────────
### MAIZE
#maize <- vroom("/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/fst/Fst_Indian_Cy0_Cy14.weir.fst")
maize <- vroom("/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/fst/new_data/Fst_Indian_Cy0_Cy14_single.weir.fst")

mai

maize <- maize %>%
  mutate(CHROM = case_when(
    CHROM == "NC_050096.1" ~ "chr1",
    CHROM == "NC_050097.1" ~ "chr2",
    CHROM == "NC_050098.1" ~ "chr3",
    CHROM == "NC_050099.1" ~ "chr4",
    CHROM == "NC_050100.1" ~ "chr5",
    CHROM == "NC_050101.1" ~ "chr6",
    CHROM == "NC_050102.1" ~ "chr7",
    CHROM == "NC_050103.1" ~ "chr8",
    CHROM == "NC_050104.1" ~ "chr9",
    CHROM == "NC_050105.1" ~ "chr10",
    TRUE ~ CHROM # Retain other values unchanged
  ))

maize <- maize %>% 
  filter(!grepl("^NW", CHROM),                  # drop scaffolds
         complete.cases(.)) %>%                # drop NAs
  mutate(WEIR_AND_COCKERHAM_FST =               # set negatives to zero
           pmax(WEIR_AND_COCKERHAM_FST, 0))

# chromosomes you want to keep
chromosomes <- paste0("chr", 1:10)

# ──────────────────────────────────────────────
# 2.  Spline windows for every chromosome ------
# ──────────────────────────────────────────────
chr_data <- lapply(chromosomes, function(chr){
  
  chr_tbl <- maize  %>% filter(CHROM == chr)  %>% arrange(POS)
  
  out <- splineAnalyze(
    Y          = chr_tbl$WEIR_AND_COCKERHAM_FST,
    map        = chr_tbl$POS,
    smoothness = 5e5,      # 500 kb
    plotRaw    = FALSE,
    plotWindows= FALSE,
    method     = 4)        # ordinary CV
  
  out$windowData |>
    mutate(CHROM = chr)
})
combined_df <- bind_rows(chr_data)


# ──────────────────────────────────────────────
# 3.  Build global coordinates   ----------------
#     (do this **after** any further filtering)
# ──────────────────────────────────────────────
plot_df <- combined_df %>% 
  filter(SNPcount > 0,                           # drop empty windows
         complete.cases(.))

# Calculate the threshold for the top 1% of MeanY
top_1_percent_threshold <- quantile(plot_df$MeanY, 0.99, na.rm = TRUE)

# Display the threshold
top_1_percent_threshold


# chromosome length & cumulative start
chrom_info <- plot_df %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(WindowStop), .groups = "drop") %>%
  mutate(cum_start = lag(cumsum(chr_len), default = 0),
         chr_mid   = cum_start + chr_len/2)

# add global positions
plot_df <- plot_df %>% 
  left_join(chrom_info, by = "CHROM") %>%
  mutate(start_g = WindowStart + cum_start,
         end_g   = WindowStop  + cum_start,
         mid_g   = (WindowStart + WindowStop)/2 + cum_start)

# ──────────────────────────────────────────────
# 4.  Rescale W‑stat to plot on FST axis  -------
# ──────────────────────────────────────────────
rng_f <- range(plot_df$MeanY,  na.rm = TRUE)
rng_w <- range(plot_df$Wstat,  na.rm = TRUE)
a <- diff(rng_f)/diff(rng_w);  b <- rng_f[1] - a*rng_w[1]
plot_df <- mutate(plot_df, Wscaled = Wstat*a + b)

# ──────────────────────────────────────────────
# 5.  Plotting theme (Nature‑ish)  --------------
# ──────────────────────────────────────────────
nature_theme <- function(base = 11){
  theme_minimal(base_size = base, base_family = "Arial") +
    theme(text            = element_text(colour = "black"),
          panel.grid      = element_blank(),
          axis.line       = element_line(colour = "black", size = .35),
          axis.ticks      = element_line(colour = "black", size = .35),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

# ──────────────────────────────────────────────
# 6A.  Upper panel  (FST line + W points) -------
# ──────────────────────────────────────────────
g_top <- ggplot(plot_df, aes(x = mid_g/1e6, group = CHROM)) +
  geom_line(aes(y = MeanY),
            colour = "#2686C1", linewidth = .45) +
  geom_point(aes(y = Wscaled),
             colour = "#E94F37", size = .9, alpha = .5) +
  scale_x_continuous(name   = "Genomic position (Mb)",
                     breaks = chrom_info$chr_mid/1e6,
                     labels = chrom_info$CHROM,
                     expand = c(0.01,0.01)) +
  scale_y_continuous(name = expression(F[ST]),
                     sec.axis = sec_axis(~ (. - b)/a, name = "W‑stat")) +
  nature_theme()

# ──────────────────────────────────────────────
# 6B.  SNP‑density strip  -----------------------
# ──────────────────────────────────────────────
g_bottom <- ggplot(plot_df) +
  geom_rect(aes(xmin = start_g/1e6, xmax = end_g/1e6,
                ymin = 0, ymax = 1, fill = SNPcount),
            colour = NA) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "SNP\ndensity",
    trans = "sqrt",
    breaks = c(min(plot_df$SNPcount), 
               median(plot_df$SNPcount), 
               max(plot_df$SNPcount)),
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth = unit(0.3, "cm"))
  ) +
  scale_x_continuous(limits = range(plot_df$end_g/1e6),
                     expand = c(0.01,0.01)) +
  nature_theme() +
  theme(axis.text   = element_blank(),
        axis.title  = element_blank(),
        legend.title = element_text(angle = 90, vjust = .5),
        plot.margin = margin(t = -5))     # squeeze against top panel

# ──────────────────────────────────────────────
# 7.  Combine & save  ---------------------------
# ──────────────────────────────────────────────
final_plot <- (g_top / g_bottom) +
  plot_layout(heights = c(4, .6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

quartz()
final_plot

ggsave("Figure_FST_Wstat_density.tiff",
       final_plot, width = 8.5, height = 5,
       dpi = 600, compression = "lzw")
