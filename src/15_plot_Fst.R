# install.packages(c("vroom","dplyr","readr","ggplot2","IRanges","GenomicRanges"), Ncp=1)
library(vroom)
library(dplyr)
library(readr)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr) 

# 1) Read all your window files
indir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window"
# Jarvis
jarvis <- list.files(indir, pattern="^chr[0-9]+_J01_J14\\.win100k_step10k\\.clean\\.tsv$", full.names = TRUE)


jarvis_annotate <- read.csv("results/Fst/Fst_top1pct_windows_with_genes_J_grain.csv")


indir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Indian_Chief/sliding_window"
# Indian Chief
ic <- list.files(indir, pattern="^chr[0-9]+_IC01_IC14\\.win100k_step10k\\.clean\\.tsv$", full.names = TRUE)

ic_annotate <- read.csv("results/Fst/Fst_top1pct_windows_with_genes_ICz_stover.csv")




# Load the files
jarvis <- vroom::vroom(jarvis, col_types = "ciiiid") %>%  # chr,start,end,midPos,Nsites,Fst
  mutate(chr = factor(chr, levels=paste0("chr", 1:10))) %>% 
  arrange(chr, start)
ic <- vroom::vroom(ic, col_types = "ciiiid") %>%  # chr,start,end,midPos,Nsites,Fst
  mutate(chr = factor(chr, levels=paste0("chr", 1:10))) %>% 
  arrange(chr, start)



glimpse(jarvis)
glimpse(ic)


glimpse(jarvis_annotate)
glimpse(ic_annotate)


# --- 1) Compute 5% cutoffs in each set ---
ic_cut   <- quantile(ic$Fst,      0.99, na.rm = TRUE)
jar_cut  <- quantile(jarvis$Fst,  0.99, na.rm = TRUE)

# (Optional) QC filter (uncomment/adjust if you want to require min sites per window)
# min_sites <- 30000
# ic     <- ic     %>% filter(Nsites >= min_sites)
# jarvis <- jarvis %>% filter(Nsites >= min_sites)

# --- 2) Keep common windows and build flags ---
both <- inner_join(
  jarvis %>% select(chr, start, end, midPos, Fst_j = Fst),
  ic      %>% select(chr, start, end, midPos, Fst_i = Fst),
  by = c("chr","start","end","midPos")
) %>%
  mutate(
    top_j = Fst_j >= jar_cut,
    top_i = Fst_i >= ic_cut,
    tag = case_when(
      top_j & top_i ~ "Both top 5%",
      top_j & !top_i ~ "Jarvis top 5% only",
      !top_j & top_i ~ "Indian Chief top 5% only",
      TRUE ~ "Neither"
    )
  )

# Quick counts + correlation (optional)
counts <- both %>% count(tag)
print(counts)
cat("Pearson r (IC vs Jarvis Fst): ",
    round(cor(both$Fst_i, both$Fst_j, use = "complete.obs"), 3), "\n")

# --- 3) Plot: Jarvis on Y, Indian Chief on X ---
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 16,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 16,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 16,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 16,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    panel.grid     = element_blank(),
    
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right","top"),
    legend.background    = element_rect(fill="white", color="grey70", size=0.4),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2,"cm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size=16),
    
    plot.margin    = margin(15, 15, 15, 15)
  )
p <- ggplot(both, aes(x = Fst_i, y = Fst_j)) +
  geom_point(aes(color = tag), size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = ic_cut,  linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = jar_cut, linetype = "dashed", linewidth = 0.4) +
  # 1:1 reference (optional)
  #geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.3) +
  scale_color_manual(
    values = c(
      "Both top 5%"              = "#D81B60",
      "Jarvis top 5% only"       = "#1E88E5",
      "Indian Chief top 5% only" = "#FFC107",
      "Neither"                  = "grey70"
    )
  ) +
  coord_equal() +
  labs(
    title = "Window FST: Indian Chief (X) vs Jarvis (Y)",
    x = "Indian Chief FST",
    y = "Jarvis FST",
    color = NULL,
    caption = paste0("Dashed lines = 95th percentile cutoffs (IC = ",
                     signif(ic_cut,3), ", Jarvis = ", signif(jar_cut,3), ")")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  plot_theme



# --- 4) Save the intersecting outliers (both top 5%) ---
# both_top5 <- both %>% filter(tag == "Both top 5%")
# write_tsv(both_top5, "FST_common_top5_both.tsv")

# (Optional) If the scatter is too dense, try a hex/2D bin view:
# ggplot(both, aes(Fst_i, Fst_j)) +
#   geom_bin2d(bins = 75) +
#   geom_vline(xintercept = ic_cut,  linetype = "dashed") +
#   geom_hline(yintercept = jar_cut, linetype = "dashed") +
#   coord_equal() + theme_bw() + labs(x="Indian Chief FST", y="Jarvis FST")






# --- your gene list (deduped) ---
genes <- unique(c(
  "Zm00001eb201190",
  "Zm00001eb201210",
  "Zm00001eb099180",
  "Zm00001eb001840",
  "Zm00001eb381160",
  "Zm00001eb387170",
  "Zm00001eb129680"
))

library(dplyr)
# library(ggrepel)

# 1) best window per gene in each pop (from earlier)
# best_J and best_I already exist from the previous code block.
# If not, rebuild them exactly as before.

# 2) all windows grid (IC ↔ Jarvis FST pairs)
both_windows <- inner_join(
  jarvis %>% select(chr, midPos, Fst_j = Fst),
  ic      %>% select(chr, midPos, Fst_i = Fst),
  by = c("chr","midPos")
)

# 3) population category per gene (from earlier)
annot_tbl <- full_join(best_J, best_I, by = "gene_id") %>%
  mutate(category = case_when(
    !is.na(J_Fst) & !is.na(IC_Fst) ~ "Common",
    !is.na(J_Fst) &  is.na(IC_Fst) ~ "Jarvis-only",
    is.na(J_Fst)  & !is.na(IC_Fst) ~ "IndianChief-only",
    TRUE ~ "NA"
  )) %>%
  select(gene_id, category)

# 4) get plot coords for EVERY gene using each pop’s best window,
#    then keep one point per gene by category
pts_J <- best_J %>%
  left_join(both_windows, by = c("J_chr"="chr", "J_midPos"="midPos")) %>%
  transmute(gene_id, source = "Jarvis-best", x = Fst_i, y = Fst_j)

pts_I <- best_I %>%
  left_join(both_windows, by = c("IC_chr"="chr", "IC_midPos"="midPos")) %>%
  transmute(gene_id, source = "IC-best", x = Fst_i, y = Fst_j)

gene_points_all <- bind_rows(pts_J, pts_I) %>%
  inner_join(annot_tbl, by = "gene_id") %>%
  mutate(rank_score = dplyr::case_when(
    category == "Common" ~ x + y,                         # pick the stronger combined point
    category == "Jarvis-only" ~ if_else(source=="Jarvis-best", 1e6, -1e6),
    category == "IndianChief-only" ~ if_else(source=="IC-best", 1e6, -1e6),
    TRUE ~ x + y
  )) %>%
  group_by(gene_id) %>%
  slice_max(order_by = rank_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# 5) add to your existing scatter `p`
p_genes <- p +
  geom_point(data = gene_points_all, aes(x = x, y = y), size = 2.5, shape = 21, stroke = 0.7, fill = NA) +
  ggrepel::geom_text_repel(data = gene_points_all,
                           aes(x = x, y = y, label = gene_id),
                           size = 3, max.overlaps = 50) +
  guides(
    color = guide_legend(
      ncol = 1,
      keyheight = unit(4, "pt"),
      keywidth  = unit(10, "pt"),
      override.aes = list(size = 2, alpha = 1)  # small legend points
    )
  ) +
  theme(
    # (A) keep base sizes modest
    text = element_text(size = 12),
    # (B) shrink legend parts
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(6, "pt"),
    legend.spacing.y  = unit(1, "pt"),
    legend.margin     = margin(2, 2, 2, 2),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    # optional: tuck it tighter into the corner
    legend.position   = c(0.95, 0.95),
    legend.justification = c("right","top")
  )

quartz()
p_genes

# Save the labeled plot
ggsave("Fst_scatter_with_genes_labeled.png", p_genes,
       width = 8, height = 6, units = "in", dpi = 300)
# 6) (optional) export a small table for the labeled genes
# ic_cut  and jar_cut are from before
# readr::write_tsv(
#   gene_points_all %>%
#     mutate(IC_Fst = x, Jarvis_Fst = y,
#            IC_top5 = IC_Fst >= ic_cut, Jarvis_top5 = Jarvis_Fst >= jar_cut) %>%
#     select(gene_id, category, IC_Fst, Jarvis_Fst, IC_top5, Jarvis_top5),
#   "genes_on_scatter.tsv"
# )
