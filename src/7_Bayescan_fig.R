################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LANDRACES MAIZE AND SORGHUM
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################


library(dplyr)
library(tidyr)

################################################################################
### Read the data
################################################################################

############# MAIZE
#### NITROGEN
bayescan_n_maize <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/maize/nitrogen/romerro_N_bayescan_result_fst.txt")

#### NHx
bayescan_n_maize <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/maize/NHx/romerro_NHx_bayescan_result_fst.txt")

#### Elevation
bayescan_n_maize <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/maize/elevation/Elevation_bayescan_result_fst.txt")

# Get the snp data
# bcftools query -f '%ID\n' all_maize2.biallelic.vcf.gz > snp_ids.txt

#####  Read the snp data
snp_data <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/maize/nitrogen/snp_ids.txt")
head(snp_data)



############# SORGHUM
#### NITROGEN
bayescan_n_maize <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/nitrogen/lasky_N_bayescan_result_fst.txt")
snp_data <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/nitrogen/snp_ids.txt")

#### NHx
bayescan_n_maize <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/NHx/lasky_NHx_bayescan_result_fst.txt")

snp_data <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/NHx/snp_ids.txt")
head(snp_data)

# Only for sorghum
drop_snp <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/nitrogen/dropped_rows.txt")
drop_snp <- read.table("/Users/nirwantandukar/Documents/Research/results/envGWAS/Bayescan/sorghum/NHx/dropped_rows.txt")

rows_to_drop <- drop_snp$V1          # numeric vector of one‑based indices

## keep everything *except* those rows
snp_data          <- snp_data[-rows_to_drop, , drop = FALSE]


################################################################################
### DATA PROCESSING
################################################################################

# Add SNP names as a new column to BayeScan results
bayescan_n_maize$snp_id <- snp_data$V1

# Separate snp_id into chr and pos
bayescan_n_maize <- bayescan_n_maize %>%
  separate(snp_id, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  mutate(chr = gsub("S", "", chr),
         pos = as.integer(pos))
head(bayescan_n_maize)

# Replace 0 qval with half of the smallest non-zero qval
min_qval <- min(bayescan_n_maize$qval[bayescan_n_maize$qval > 0], na.rm = TRUE)
bayescan_n_maize$qval[bayescan_n_maize$qval == 0] <- min_qval / 2

# find the largest finite log10.PO. (<1000) and the smallest non‑zero q‑value
max_PO  <- max(bayescan_n_maize$log10.PO.[bayescan_n_maize$log10.PO. < 1000], na.rm = TRUE)
min_q   <- min(bayescan_n_maize$qval[bayescan_n_maize$qval > 0],                na.rm = TRUE)

# replace 1000 with 2×max_PO, and 0 with ½×min_q
bayescan_n_maize <- transform(
  bayescan_n_maize,
  log10.PO. = replace(log10.PO., log10.PO. == 1000, max_PO),
  qval      = replace(qval,      qval      ==    0, min_q  / 2)
)


# Dummy variable for plotting
plot_bayesian <- bayescan_n_maize

# Filter based on pvalue jeffreys scale
#bayescan_n_maize <- subset(bayescan_n_maize, qval < 0.05 & prob > 0.76)
head(bayescan_n_maize)


str(plot_bayesian)

################################################################################
### PLOTTING
################################################################################

# Ensure all chromosomes are present in plot_bayesian
plot_bayesian <- plot_bayesian %>%
  mutate(chr = factor(chr, levels = as.character(1:10)))

# Compute chromosome offsets (modified to handle missing chromosomes)
chr_info <- plot_bayesian %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

# Handle cases where some chromosomes might have no data
if(nrow(chr_info) < 10) {
  missing_chrs <- setdiff(as.character(1:10), chr_info$chr)
  chr_info <- bind_rows(
    chr_info,
    tibble(chr = missing_chrs, chr_len = NA, offset = NA)
  ) %>% arrange(chr)
}

# Join offset
plot_df <- plot_bayesian %>%
  left_join(chr_info, by = "chr") %>%
  mutate(pos_cum = pos + offset)

# Create axis labels only for chromosomes that have data
axis_df <- plot_df %>%
  group_by(chr) %>%
  summarise(center = mean(range(pos_cum, na.rm = TRUE)), .groups = "drop") %>%
  filter(!is.na(center))


# Calculate the top 1% FST threshold
top_1_percent_threshold <- quantile(plot_bayesian$fst, probs = 0.99, na.rm = TRUE)

# Print the threshold value
print(paste("Top 1% FST threshold:", top_1_percent_threshold))

# Manhattan plot
manhattan_plot <- ggplot(plot_df, aes(x = pos_cum, y = fst)) +
  
  # Use geom_jitter instead of geom_point for better performance with many points
  geom_jitter(aes(color = as.factor(as.numeric(as.character(chr)) %% 2)), 
              alpha = 0.6, 
              size = 0.8, 
              width = 0.1, 
              height = 0) +
  
  # Simplified color scale
  scale_color_manual(values = c("#1F78B4", "#A6CEE3")) +
  
  # Threshold line (simplified)
  geom_hline(yintercept = top_1_percent_threshold, 
             linetype = "dashed", 
             color = "#E31A1C",
             linewidth = 0.5) +
  
  # Lightweight text annotation instead of geom_label
  annotate("text",
           x = max(plot_df$pos_cum, na.rm = TRUE) * 0.95,
           y = top_1_percent_threshold * 1.05,
           label = paste("Top 1% =", round(top_1_percent_threshold, 3)),
           color = "#E31A1C",
           size = 3) +
  
  # Chromosome axis (simplified)
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chr,
    expand = c(0.01, 0.01)
  ) +
  
  # Lightweight theme
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(F[ST]),
    title = "Genome-wide FST Manhattan Plot"
  )

# Display with optimal dimensions
quartz(width = 12, height = 6)
print(manhattan_plot)

# Save as high-quality PDF
ggsave("FST_Manhattan_plot_elevation_maize.pdf", 
       plot = manhattan_plot, 
       width = 12, 
       height = 6, 
       device = cairo_pdf)




# Plot Fst vs qvalue
#plot_bayesian$alpha_sign <- ifelse(plot_bayesian$alpha > 0, "diversifying",
#                                   ifelse(plot_bayesian$alpha < 0, "balancing", "neutral"))
plot_bayesian$alpha_sign <- ifelse(plot_bayesian$alpha > 0, "diversifying", "balancing")
                                   


# First, ensure alpha_sign is properly factored (in case you want neutral later)
# plot_bayesian <- plot_bayesian %>%
#   mutate(alpha_sign = factor(alpha_sign, 
#                              levels = c("diversifying", "balancing", "neutral"),
#                              labels = c("Diversifying", "Balancing", "Neutral")))

plot_bayesian <- plot_bayesian %>%
  mutate(alpha_sign = factor(alpha_sign,
                             levels = c("diversifying", "balancing"),
                             labels = c("Diversifying", "Balancing")))


# Create the plot with publication-quality styling
library(ggplot2)
selection_plot <- ggplot(plot_bayesian, aes(x = log10.PO., y = fst, color = alpha_sign)) +
  geom_point(alpha = 0.8, size = 2, shape = 16) +  # Adjust point aesthetics
  scale_color_manual(
    name = "Selection Type",
    values = c(
      "Diversifying" = "#2E8B57",  # More professional green (SeaGreen)
      "Balancing" = "#D55E00",     # Professional orange
      "Neutral" = "gray70"         # Neutral color
    ),
    drop = FALSE  # Show all categories even if not present
  ) +
  # Add significance threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", linewidth = 0.5) +
  # Add top 1% FST threshold if desired
  geom_hline(yintercept = quantile(plot_bayesian$fst, 0.99, na.rm = TRUE), 
             linetype = "dotted", color = "black", linewidth = 0.5) +
  # Labels and titles
  labs(
    x = expression(log[10]("Posterior Odds")),
    y = expression(F[ST]),
    title = "Selection Signature Analysis",
    subtitle = "Bayescan results for soil nitrogen in maize landraces colored by selection type"
  ) +
  # Publication-quality theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  # Adjust scales if needed
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))

# Display the plot
quartz(width = 8, height = 6)  # Optimal size for single-column figures
print(selection_plot)

# Save as high-quality PDF for publication
ggsave("Bayescan_selection_plot_Nitrogen_maize.png", 
       plot = selection_plot, 
       width = 8, 
       height = 6, 
       bg = "white",
       device = png)




# Load the GFF annotation file
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")

# Filter for genes only
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Sorghum
# Load your reference GFF file - replace with the actual path
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]




# Convert to GRanges
bayescan_gr <- GRanges(
  # In maize
  seqnames = paste0("chr", bayescan_n_maize$chr),
  # In sorghum no chr
  #seqnames = paste0(bayescan_n_maize$chr),
  ranges = IRanges(start = bayescan_n_maize$pos, end = bayescan_n_maize$pos),
  prob = bayescan_n_maize$prob,
  qval = bayescan_n_maize$qval,
  alpha = bayescan_n_maize$alpha,
  fst = bayescan_n_maize$fst
)

# Find overlapping or nearby genes within 50,000 bp of high Fst SNPs
overlapping_genes <- findOverlaps(bayescan_gr, genes_only, maxgap = 25000)


# Extract gene information
genes_near_high_fst <- genes_only[subjectHits(overlapping_genes)]
gene_names <- mcols(genes_near_high_fst)$ID  # Replace 'ID' with the appropriate column for gene IDs

unique(gene_names)

# Convert GRanges to data frame
genes_df <- as.data.frame(genes_near_high_fst)

# Extract SNP info that overlaps with genes (query hits)
snp_hits <- bayescan_n_maize[queryHits(overlapping_genes), ]

# Combine SNP and gene info
combined_info <- cbind(
  chr         = snp_hits$chr,
  snp_pos     = snp_hits$pos,
  snp_id      = snp_hits$snp_id,
  qval        = snp_hits$qval,
  alpha       = snp_hits$alpha,
  fst         = snp_hits$fst,
  log10_PO    = snp_hits$log10.PO.,
  gene_chr    = as.character(genes_df$seqnames),
  gene_start  = genes_df$start,
  gene_end    = genes_df$end,
  strand      = as.character(genes_df$strand),
  #gene_id     = genes_df$gene_id,
  gene_name   =  genes_df$ID,
  description = genes_df$description,
  biotype     = genes_df$biotype
)



# Sorghum
combined_info <- cbind(
  chr         = snp_hits$chr,
  snp_pos     = snp_hits$pos,
  snp_id      = snp_hits$snp_id,
  qval        = snp_hits$qval,
  alpha       = snp_hits$alpha,
  fst         = snp_hits$fst,
  log10_PO    = snp_hits$log10.PO.,
  gene_chr    = as.character(genes_df$seqnames),
  gene_start  = genes_df$start,
  gene_end    = genes_df$end,
  strand      = as.character(genes_df$strand),
  gene_id     = genes_df$gene_id,
  gene_name   = sub("gene:", "", genes_df$ID),
  description = genes_df$description,
  biotype     = genes_df$biotype
)



# Optional: remove "gene:" prefix from gene_name
# Maize
combined_info <- combined_info[,c(12,1:7)]

# Sorghum
combined_info <- combined_info[,c(13,1:7)]
head(combined_info)


# Remove duplicates by gene_name
combined_info <- as.data.frame(combined_info)
combined_info_unique <- combined_info[!duplicated(combined_info$gene_name), ]

# Qval < 0.05
combined_info <- combined_info[combined_info$qval < 0.05, ]

# Save to file
write.csv(combined_info, "Bayescan_Elevation_romerror_maize_25kextend.csv", row.names = FALSE , quote = FALSE)
#write.csv(combined_info, "Bayescan_NHx_lasky_sorghum_25kextend.csv", row.names = FALSE , quote = FALSE)
getwd()




#### Plotting ####
library(dplyr)
library(ggplot2)

# Step 1: Prepare cumulative position
plot_bayesian <- plot_bayesian %>%
  arrange(chr, pos) %>%
  group_by(chr) %>%
  mutate(pos_rel = pos) %>%
  ungroup()

# Step 2: Calculate cumulative positions
chr_offsets <- plot_bayesian %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(chr_len = max(pos)) %>%
  dplyr::mutate(offset = lag(cumsum(chr_len), default = 0))

plot_bayesian <- plot_bayesian %>%
  dplyr::left_join(chr_offsets, by = "chr") %>%
  dplyr::mutate(pos_cum = pos + offset,
                log_qval = -log10(qval))

# Step 3: Axis labels in the center of each chromosome
axis_df <- plot_bayesian %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(center = (min(pos_cum) + max(pos_cum)) / 2)

plot_bayesian_clean <- plot_bayesian %>% 
  filter(is.finite(log_qval))

axis_df <- plot_bayesian_clean %>%
  dplyr::group_by(chr) %>%
  dplyr::summarize(center = mean(pos_cum, na.rm = TRUE))


# Step 4: Plot
quartz()
ggplot(plot_bayesian_clean, aes(x = pos_cum, y = log_qval, color = as.factor(as.numeric(chr) %% 2))) +
  geom_point(size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "gray50")) +
  scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
  labs(x = "Chromosome", y = expression(-log[10](qvalue)), title = "Bayescan Manhattan Plot") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )




# Fst vs qvalue


# Your data frame: plot_bayesian
library(ggplot2)

# Calculate -log10(qvalue) if not already present
plot_bayesian$log_qval <- -log10(plot_bayesian$qval)

# Plot
quartz()
ggplot(plot_bayesian, aes(x = log_qval, y = fst)) +
  geom_point(aes(color = qval < 0.05), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(
    x = expression(-log[10](qvalue)),
    y = expression(F[ST]),
    title = "BayeScan FST vs -log10(q-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



## --- classify BayeScan SNPs and make “balancing / diversifying / neutral” plot ----
library(dplyr)
library(ggplot2)

## 1. define thresholds -----------------------------------------------------------
fdr_cutoff <- 0.05      # q‑value threshold (FDR)
po_cutoff  <- 2         # log10(PO) threshold for the reference line

## 2. cap extreme log10(PO) values (Bayescan uses 1000 for ∞) ---------------------
max_PO <- max(plot_bayesian$log10.PO.[plot_bayesian$log10.PO. < 1000], na.rm = TRUE)
plot_bayesian$log10.PO.[plot_bayesian$log10.PO. == 1000] <- max_PO * 1.05  # small offset

## 3. add selection category ------------------------------------------------------
quartz()
plot_bayesian <- plot_bayesian %>%
  mutate(
    Selection = case_when(
      qval < fdr_cutoff & alpha  > 0 ~ "diversifying",
      qval < fdr_cutoff & alpha <= 0 ~ "balancing",
      TRUE                           ~ "neutral"
    ),
    Selection = factor(Selection, levels = c("balancing", "diversifying", "neutral"))
  )

## 4.plot ------------------------------------------------------------------------
quartz()            # comment out if not on macOS / using RStudio graphics device
ggplot(plot_bayesian, aes(x = log10.PO., y = fst, color = Selection)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(
    values = c(
      balancing     = "#F8766D",  # reddish
      diversifying = "#00BA38",  # green
      neutral      = "#619CFF"   # blue
    )
  ) +
  geom_vline(xintercept = 2, linetype = "solid", color = "black") +
  labs(
    x = expression(log[10](PO)),
    y = expression(F[ST]),
    color = "Selection",
    title = "BayeScan Outlier Classification"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


table(plot_bayesian$Selection)
summary(plot_bayesian$alpha[plot_bayesian$qval < 0.05])



plot_bayesian$alpha_sign <- ifelse(plot_bayesian$alpha > 0, "positive", "negative")

quartz()
ggplot(plot_bayesian, aes(x = log10.PO., y = fst, color = alpha_sign)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(positive = "#00BA38", negative = "#F8766D")) +
  geom_vline(xintercept = 2, linetype = "solid", color = "black") +
  labs(
    x = expression(log[10](PO)), y = expression(F[ST]), color = "Alpha sign",
    title = "BayeScan: Divergent vs Balancing Potential Signals"
  ) +
  theme_minimal()



















################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# INDIAN CHIEF AND JARVIS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################

################################################################################
### Read the data
################################################################################

############# Jarvis

bayescan_jarvis <- read.table("/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/bayescan/Jarvis/Jarvis_bayescan_result_fst.txt")


# Get the snp data
# This is different than the maize one
# bcftools query -f '%CHROM\t%POS\n' Indian_cycle0_filtered_0.1MAF_onlyCHR.vcf.gz > snp_ids.txt

#####  Read the snp data
snp_data <- read.table("/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/bayescan/Jarvis/snp_ids.txt")
head(snp_data)
colnames(snp_data) <- c("snp_id", "pos")
snp_data$snp <- paste0(snp_data$snp_id, "_", snp_data$pos)


# Only for sorghum
drop_snp <- read.table("/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/bayescan/Jarvis/dropped_rows_Jarvis.txt")

rows_to_drop <- drop_snp$V1          # numeric vector of one‑based indices

## keep everything *except* those rows
snp_data <- snp_data[-rows_to_drop, , drop = FALSE]



################################################################################
### DATA PROCESSING
################################################################################

# Add SNP names as a new column to BayeScan results
bayescan_jarvis$chr <- snp_data$snp_id
bayescan_jarvis$pos <- snp_data$pos

# Replace 0 qval with half of the smallest non-zero qval
min_qval <- min(bayescan_jarvis$qval[bayescan_jarvis$qval > 0], na.rm = TRUE)
bayescan_jarvis$qval[bayescan_jarvis$qval == 0] <- min_qval / 2

# find the largest finite log10.PO. (<1000) and the smallest non‑zero q‑value
max_PO  <- max(bayescan_jarvis$log10.PO.[bayescan_jarvis$log10.PO. < 1000], na.rm = TRUE)
min_q   <- min(bayescan_jarvis$qval[bayescan_jarvis$qval > 0],                na.rm = TRUE)

# replace 1000 with 2×max_PO, and 0 with ½×min_q
bayescan_jarvis <- transform(
  bayescan_jarvis,
  log10.PO. = replace(log10.PO., log10.PO. == 1000, max_PO),
  qval      = replace(qval,      qval      ==    0, min_q  / 2)
)


# Dummy variable for plotting
plot_bayesian <- bayescan_jarvis

# Filter based on pvalue jeffreys scale
#bayescan_n_maize <- subset(bayescan_n_maize, qval < 0.05 & prob > 0.76)
head(bayescan_jarvis)


str(plot_bayesian)




################################################################################
### PLOTTING
################################################################################

# Ensure all chromosomes are present in plot_bayesian
plot_bayesian <- plot_bayesian %>%
  mutate(chr = factor(chr, levels = as.character(1:10)))

# Compute chromosome offsets (modified to handle missing chromosomes)
chr_info <- plot_bayesian %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

# Handle cases where some chromosomes might have no data
if(nrow(chr_info) < 10) {
  missing_chrs <- setdiff(as.character(1:10), chr_info$chr)
  chr_info <- bind_rows(
    chr_info,
    tibble(chr = missing_chrs, chr_len = NA, offset = NA)
  ) %>% arrange(chr)
}

# Join offset
plot_df <- plot_bayesian %>%
  left_join(chr_info, by = "chr") %>%
  mutate(pos_cum = pos + offset)

# Create axis labels only for chromosomes that have data
axis_df <- plot_df %>%
  group_by(chr) %>%
  summarise(center = mean(range(pos_cum, na.rm = TRUE)), .groups = "drop") %>%
  filter(!is.na(center))


# Calculate the top 1% FST threshold
top_1_percent_threshold <- quantile(plot_bayesian$fst, probs = 0.99, na.rm = TRUE)

# Print the threshold value
print(paste("Top 1% FST threshold:", top_1_percent_threshold))

# Manhattan plot
manhattan_plot <- ggplot(plot_df, aes(x = pos_cum, y = fst)) +
  
  # Use geom_jitter instead of geom_point for better performance with many points
  geom_jitter(aes(color = as.factor(as.numeric(as.character(chr)) %% 2)), 
              alpha = 0.6, 
              size = 0.8, 
              width = 0.1, 
              height = 0) +
  
  # Simplified color scale
  scale_color_manual(values = c("#1F78B4", "#A6CEE3")) +
  
  # Threshold line (simplified)
  geom_hline(yintercept = top_1_percent_threshold, 
             linetype = "dashed", 
             color = "#E31A1C",
             linewidth = 0.5) +
  
  # Lightweight text annotation instead of geom_label
  annotate("text",
           x = max(plot_df$pos_cum, na.rm = TRUE) * 0.95,
           y = top_1_percent_threshold * 1.05,
           label = paste("Top 1% =", round(top_1_percent_threshold, 3)),
           color = "#E31A1C",
           size = 3) +
  
  # Chromosome axis (simplified)
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chr,
    expand = c(0.01, 0.01)
  ) +
  
  # Lightweight theme
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(F[ST]),
    title = "Genome-wide FST Manhattan Plot"
  )

# Display with optimal dimensions
quartz(width = 12, height = 6)
print(manhattan_plot)

# Save as high-quality PDF
ggsave("FST_Manhattan_plot_Jarvis.pdf", 
       plot = manhattan_plot, 
       width = 12, 
       height = 6, 
       device = cairo_pdf)




# Plot Fst vs qvalue
plot_bayesian$alpha_sign <- ifelse(plot_bayesian$alpha > 0, "diversifying",
                                   ifelse(plot_bayesian$alpha < 0, "balancing", "neutral"))

# First, ensure alpha_sign is properly factored (in case you want neutral later)
plot_bayesian <- plot_bayesian %>%
  mutate(alpha_sign = factor(alpha_sign, 
                             levels = c("diversifying", "balancing", "neutral"),
                             labels = c("Diversifying", "Balancing", "Neutral")))

# Create the plot with publication-quality styling
library(ggplot2)
selection_plot <- ggplot(plot_bayesian, aes(x = log10.PO., y = fst, color = alpha_sign)) +
  geom_point(alpha = 0.8, size = 2, shape = 16) +  # Adjust point aesthetics
  scale_color_manual(
    name = "Selection Type",
    values = c(
      "Diversifying" = "#2E8B57",  # More professional green (SeaGreen)
      "Balancing" = "#D55E00",     # Professional orange
      "Neutral" = "gray70"         # Neutral color
    ),
    drop = FALSE  # Show all categories even if not present
  ) +
  # Add significance threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", linewidth = 0.5) +
  # Add top 1% FST threshold if desired
  geom_hline(yintercept = quantile(plot_bayesian$fst, 0.99, na.rm = TRUE), 
             linetype = "dotted", color = "black", linewidth = 0.5) +
  # Labels and titles
  labs(
    x = expression(log[10]("Posterior Odds")),
    y = expression(F[ST]),
    title = "Selection Signature Analysis",
    subtitle = "Bayescan results colored by selection type"
  ) +
  # Publication-quality theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  # Adjust scales if needed
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))

# Display the plot
quartz(width = 8, height = 6)  # Optimal size for single-column figures
print(selection_plot)

# Save as high-quality PDF for publication
ggsave("Bayescan_selection_plot_NHx_sorghum.pdf", 
       plot = selection_plot, 
       width = 8, 
       height = 6, 
       device = cairo_pdf)

