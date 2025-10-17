#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
# PCA analysis
################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(SNPRelate)

###PCA analysis
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting
setwd("/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/")

# Converting vcf to gds
vcf_ij <- snpgdsVCF2GDS("Indian_Jarvis_chr10_renamed.vcf.gz", "Indian_Jarvis_chr10.gds", method = "copy.num.of.ref")

# get the gds fil
gdsfile <- snpgdsOpen("Indian_Jarvis_chr10.gds")


##LD-based SNP pruning
set.seed(1000)

# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
snpset <- snpgdsLDpruning(gdsfile, ld.threshold = 0.2)

# Get all the SNPs
snpset.id <- unlist(unname(snpset))

# Run PCA
pca <- snpgdsPCA(gdsfile, snp.id = snpset.id, num.thread = 10)

# Get the eigenvalues
pca$eigval[1]

library(ggplot2)

# Convert PCA results to a data frame
pc_df <- data.frame(
  sample.id = pca$sample.id,  # Sample IDs
  PC1 = pca$eigenvect[,1],    # First Principal Component
  PC2 = pca$eigenvect[,2]     # Second Principal Component
)

# View the first few rows
head(pc_df)

# Add a column to categorize samples by prefix
pc_df$group <- case_when(
  grepl("^I01", pc_df$sample.id) ~ "I0",
  grepl("^I14", pc_df$sample.id) ~ "I14",
  grepl("^J01", pc_df$sample.id) ~ "J0",
  grepl("^J14", pc_df$sample.id) ~ "J14",
  TRUE ~ "Other"
)


# PCA Scatter Plot
# PCA Plot with Labels by Group
# PCA Plot without dot labels and fixed legend
pca_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8, size = 3, stroke = 0.5) +  # Slightly larger points with subtle stroke
  theme_minimal(base_size = 18) +
  labs(
    x = paste0("PC 1 (", round(pca$varprop[1] * 100, 2), "% variance)"),
    y = paste0("PC 2 (", round(pca$varprop[2] * 100, 2), "% variance)"),
    title = "PCA of Indian Chief and Jarvis Lines",
    subtitle = "Comparison between generations 0 and 14"
  ) +
  scale_color_manual(
    name = "Lines and Generation",
    values = c("I0" = "#E41A1C", "I14" = "#377EB8", "J0" = "#4DAF4A", "J14" = "#984EA3"),
    labels = c("I0" = "Indian Chief G0", "I14" = "Indian Chief G14", 
               "J0" = "Jarvis G0", "J14" = "Jarvis G14")
  ) +
  theme(
    text = element_text(family = "Arial"),  # Use a clean font
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +  # Make legend symbols more visible
  coord_fixed(ratio = 1)  # Keep aspect ratio square for proper PCA interpretation

quartz()
pca_plot



pca_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8, size = 3, stroke = 0.5) +
  theme_minimal(base_size = 18) +
  labs(
    x = paste0("PC 1 (", round(pca$varprop[1] * 100, 2), "% variance)"),
    y = paste0("PC 2 (", round(pca$varprop[2] * 100, 2), "% variance)"),
    title = "PCA of Indian Chief and Jarvis Lines",
    subtitle = "Comparison between generations 0 and 14"
  ) +
  scale_color_manual(
    name = "Line and Generation",
    values = c("I0" = "#E41A1C", "I14" = "#377EB8", "J0" = "#4DAF4A", "J14" = "#984EA3"),
    labels = c("I0" = "Indian Chief G0", "I14" = "Indian Chief G14", 
               "J0" = "Jarvis G0", "J14" = "Jarvis G14")
  ) +
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),  # Slightly smaller for space
    legend.title = element_text(size = 14, face = "bold", color = "black"),  # Slightly smaller
    legend.position = c(0.95, 0.05),  # Bottom right position
    legend.justification = c(1, 0),  # Anchored to bottom right
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key = element_blank(),
    legend.spacing = unit(0.2, "cm"),
    legend.box.margin = margin(5, 5, 5, 5),  # Small margin around legend
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  # Slightly smaller legend symbols
  coord_fixed(ratio = 1)


# Save the plot
ggsave("PCA_Indian_Jarvis.png", plot = pca_plot, width = 10, height = 8, dpi = 300)

