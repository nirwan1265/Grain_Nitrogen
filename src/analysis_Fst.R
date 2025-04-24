################################################################################
###### LOAD PACKAGES
################################################################################

# install.packages(c("vcfR","vegan","ggplot2","reshape2"))
library(vcfR); library(vegan); library(ggplot2); library(reshape2); library(vroom); library(dplyr); library(GenWin)


################################################################################
###### DATA DIVISION
################################################################################

# Get the phenotype
maize_n <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_nitrogen.csv") %>% dplyr::select(c(1,2))
maize_n <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_inorg_N.csv") %>% dplyr::select(c(1,3))
maize_n <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_OlsenP_all_data.csv") %>% dplyr::select(c(1,2))


maize_n <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/sorghum_nitrogen.csv") %>% dplyr::select(c(1,2))
maize_n <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/sorghum_OlsenP_TP.csv") %>% dplyr::select(c(1,3))
maize_n <- maize_n[complete.cases(maize_n), ] # remove rows with NA

# Sort and compute quantiles
phen_col <- maize_n[,2]
q1 <- quantile(phen_col, 0.25, na.rm = TRUE)
q3 <- quantile(phen_col, 0.75, na.rm = TRUE)

# Split into two groups
low_group  <- maize_n %>% filter(.[[2]] <= q1)
high_group <- maize_n %>% filter(.[[2]] >= q3)


################################################################################
###### SUBSET QUARTILES
################################################################################

# Set working directory to where the VCF is
vcf_dir <- "/Users/nirwantandukar/Documents/Research/data/maize_genotype/RomeroNavarro"
vcf_file <- file.path(vcf_dir, "all_maize2.vcf.gz")

vcf_dir <- "/Users/nirwantandukar/Documents/Research/data/sorghum_genotype/Lasky/"
vcf_file <- file.path(vcf_dir, "sorghum_lasky_allchrom2.vcf.gz")



# Create sample lists for low and high N groups
writeLines(low_group[[1]],  file.path(vcf_dir, "low_N_samples.txt"))
writeLines(high_group[[1]], file.path(vcf_dir, "high_N_samples.txt"))

# Create filtered VCFs using bcftools
low_out  <- file.path(vcf_dir, "low_N_3rd_quartile.vcf.gz")
high_out <- file.path(vcf_dir, "high_N_1st_quartile.vcf.gz")

# Run bcftools from R
system(paste("bcftools view -Oz -S", 
             shQuote(file.path(vcf_dir, "low_N_samples.txt")), 
             "-o", shQuote(low_out), 
             shQuote(vcf_file)))

system(paste("bcftools view -Oz -S", 
             shQuote(file.path(vcf_dir, "high_N_samples.txt")), 
             "-o", shQuote(high_out), 
             shQuote(vcf_file)))

# Index the output files
system(paste("bcftools index", shQuote(low_out)))
system(paste("bcftools index", shQuote(high_out)))


# Pick chromosome (change to whatever you want)
chr <- "10"

# Define new output paths
low_chr1  <- file.path(vcf_dir, paste0("low_N_chr", chr, ".vcf.gz"))
high_chr1 <- file.path(vcf_dir, paste0("high_N_chr", chr, ".vcf.gz"))

# Extract chr1
system(paste("bcftools view -Oz -r", chr, "-o", shQuote(low_chr1), shQuote(low_out)))
system(paste("bcftools view -Oz -r", chr, "-o", shQuote(high_chr1), shQuote(high_out)))

# Index new files
system(paste("bcftools index", shQuote(low_chr1)))
system(paste("bcftools index", shQuote(high_chr1)))



################################################################################
###### JACCARD SIMILARITY
################################################################################

library(vegan)
library(vcfR)

### ---- SETTINGS ----
#vcf_dir   <- "/Users/nirwantandukar/Documents/Research/data/maize_genotype/RomeroNavarro"
# Do high and low
vcf_in    <- file.path(vcf_dir, "high_N_chr10.vcf.gz")   # or low_N_3rd_quartile.vcf.gz
vcf_out   <- file.path(vcf_dir, "high_N_chr10_pruned.vcf.gz")
#cutoff    <- 0.40   # Jaccard similarity threshold
cutoff <- quantile(sim.mat[lower.tri(sim.mat)], 0.90)   # ≈ 0.34


vcf_in    <- file.path(vcf_dir, "low_N_chr10.vcf.gz")   # or low_N_3rd_quartile.vcf.gz
vcf_out   <- file.path(vcf_dir, "low_N_chr10_pruned.vcf.gz")
#cutoff    <- 0.40   # Jaccard similarity threshold
cutoff <- quantile(sim.mat[lower.tri(sim.mat)], 0.90)   # ≈ 0.34



### ---- 1. LOAD VCF AND BINARISE ----
vcf <- read.vcfR(vcf_in, verbose = FALSE)
gt  <- extract.gt(vcf)                   # strings like 0/0 0/1 1/1 ./.
bin <- apply(gt, 2, function(x) as.numeric(x != "0/0"))
bin[is.na(bin)] <- 0  

# Remove missing
#bin <- bin[!is.na(bin)]                  # remove missing
dim(gt_low)
ncol(gt_low)

### ---- 2. PAIRWISE JACCARD ----
dist.mat  <- vegdist(t(bin), method = "jaccard")
sim.mat   <- 1 - as.matrix(dist.mat)     # similarity

### ---- 3. ITERATIVE PRUNING ----
keep <- rep(TRUE, nrow(sim.mat))
while (TRUE) {
  m <- sim.mat[keep, keep]
  diag(m) <- 0
  if (max(m) < cutoff) break
  # pick the pair with highest similarity
  idx  <- which(m == max(m), arr.ind = TRUE)[1, ]
  # decide which one to drop: remove the one more similar to everyone (higher mean)
  means <- rowMeans(m)
  drop  <- if (means[idx[1]] >= means[idx[2]]) idx[1] else idx[2]
  keep[which(keep)[drop]] <- FALSE
}
samples_keep <- rownames(sim.mat)[keep]

# what's the maximum similarity in the matrix?
max_sim <- max(sim.mat[lower.tri(sim.mat)])
summary_sim <- summary(sim.mat[lower.tri(sim.mat)])
cat("Max pairwise Jaccard:", round(max_sim,3), "\n")
print(summary_sim)

cat(sprintf("Keeping %d of %d samples (cut‑off = %.2f)\n",
            length(samples_keep), nrow(sim.mat), cutoff))


### ---- 4. WRITE KEEP LIST & FILTER VCF ----
keep_file <- file.path(vcf_dir, "keep_highN_chr1.txt")
writeLines(samples_keep, keep_file)

keep_file <- file.path(vcf_dir, "keep_lowN_chr1.txt")
writeLines(samples_keep, keep_file)


cmd <- sprintf("bcftools view -Oz -S %s -o %s %s",
               shQuote(keep_file), shQuote(vcf_out), shQuote(vcf_in))
system(cmd)
system(sprintf("bcftools index %s", shQuote(vcf_out)))



################################################################################
###### FST CALCULATION
################################################################################

# -------- FILE LOCATIONS ------------------------------------------------------
#vcf_dir <- "/Users/nirwantandukar/Documents/Research/data/maize_genotype/RomeroNavarro"

# Full‑genome VCF
full_vcf <- file.path(vcf_dir, "all_maize2.vcf.gz")

full_vcf <- file.path(vcf_dir, "sorghum_lasky_allchrom2.vcf.gz")

# Chr10‑pruned VCFs (only used to get sample names)
low_vcf_chr10  <- file.path(vcf_dir, "low_N_chr10_pruned.vcf.gz")
high_vcf_chr10 <- file.path(vcf_dir, "high_N_chr10_pruned.vcf.gz")

# Output prefix for the FST run
out_prefix <- file.path(vcf_dir, "genomewide_low_vs_high")
fst_file   <- paste0(out_prefix, ".weir.fst")


### 1. Extract pruned sample names from the chr10 VCFs
low_pruned  <- system(sprintf("bcftools query -l %s",  shQuote(low_vcf_chr10)),  intern = TRUE)
high_pruned <- system(sprintf("bcftools query -l %s",  shQuote(high_vcf_chr10)), intern = TRUE)

# write each population’s list
low_samples_file  <- file.path(vcf_dir, "low_samples_chr10_list.txt")
high_samples_file <- file.path(vcf_dir, "high_samples_chr10_list.txt")
writeLines(low_pruned,  low_samples_file)
writeLines(high_pruned, high_samples_file)


### 2. Subset the full‑genome VCF to only those pruned samples
# combine both lists into one “all_pruned” file
all_pruned       <- unique(c(low_pruned, high_pruned))
all_pruned_file  <- file.path(vcf_dir, "all_pruned_samples.txt")
writeLines(all_pruned, all_pruned_file)

# subset full VCF
pruned_full_vcf <- file.path(vcf_dir, "pruned_all_maize2.vcf.gz")
system(sprintf(
  "bcftools view -Oz -S %s -o %s %s",
  shQuote(all_pruned_file), shQuote(pruned_full_vcf), shQuote(full_vcf)
))
system(sprintf("bcftools index %s", shQuote(pruned_full_vcf)))


### 3. Run vcftools per‑SNP FST on the pruned full‑genome VCF
cmd_fst <- sprintf(
  "vcftools --gzvcf %s \\
             --weir-fst-pop %s \\
             --weir-fst-pop %s \\
             --out %s",
  shQuote(pruned_full_vcf),
  shQuote(low_samples_file),
  shQuote(high_samples_file),
  shQuote(out_prefix)
)
system(cmd_fst)


### 4. Load & inspect the FST output
fst_df <- read.table(fst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Chromosomes in FST result:\n")
print(unique(fst_df$CHROM))
summary(fst_df$WEIR_AND_COCKERHAM_FST)


### 5. (Optional) Quick Manhattan‑style plot
quartz()
ggplot(fst_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(alpha = 0.5, size = 0.6) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 5) +
  labs(title = "Genome‑wide per‑SNP FST (Low vs High N)",
       x = "Position", y = "FST") +
  theme_minimal()






### SPLINE

# remove negative values
fst_df <- fst_df[!fst_df$WEIR_AND_COCKERHAM_FST < 0,]

chromosomes <- c(1:10)
# Initialize an empty list to store data frames
chr_data <- list()

# Loop through each chromosome
for (chr in chromosomes) {
  # Filter data for the chromosome
  chr_maize <- fst_df[which(fst_df$CHROM == chr), ]
  
  # Perform spline analysis
  chr_spline <- splineAnalyze(
    Y = chr_maize$WEIR_AND_COCKERHAM_FST,
    map = chr_maize$POS,
    smoothness = 100000,
    plotRaw = FALSE,
    plotWindows = TRUE,
    method = 4
  )
  
  # Extract window data
  chr_window <- chr_spline[["windowData"]]
  chr_window$CHROM <- chr  # Add chromosome column
  
  # Append to the list
  chr_data[[chr]] <- chr_window
}

# Combine all chromosome data into one data frame
combined_df <- do.call(rbind, chr_data)
filtered_df <- combined_df[combined_df$SNPcount >= 0, ]
filtered_df <- filtered_df[complete.cases(filtered_df), ]  # Remove rows with NA


# Ensure chromosomes are sorted in proper order
#filtered_df$CHROM <- factor(filtered_df$CHROM, levels = paste0("chr", c(1:10)))

# Create cumulative positions for each chromosome
filtered_df <- filtered_df %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(chrom_offset = lag(cumsum(WindowStop - WindowStart), default = 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cumulative_pos = WindowStart + ifelse(is.na(chrom_offset), 0, chrom_offset))

# Adjust cumulative position by adding offsets for each chromosome
chrom_offsets <- filtered_df %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarize(offset = max(cumulative_pos, na.rm = TRUE)) %>%
  dplyr::mutate(cumulative_offset = lag(cumsum(offset), default = 0))

filtered_df <- filtered_df %>%
  dplyr::left_join(chrom_offsets, by = "CHROM") %>%
  dplyr::mutate(cumulative_pos = cumulative_pos + cumulative_offset)



# Calculate the threshold for the top 1% of MeanY
top_1_percent_threshold <- quantile(filtered_df$MeanY, 0.99, na.rm = TRUE)

# Calculate the threshold for the top 5% of MeanY
top_1_percent_threshold <- quantile(filtered_df$MeanY, 0.95, na.rm = TRUE)


# Display the threshold
top_1_percent_threshold

# Manhattan plot with alternating black and gray colors
quartz()
plot_fst <- ggplot(filtered_df, aes(x = cumulative_pos, y = MeanY, color = as.factor(CHROM))) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = rep(c("black", "darkgray"), length.out = length(unique(filtered_df$CHROM)))) +
  geom_hline(yintercept = top_1_percent_threshold, linetype = "dotted", color = "red", size = 0.8) + # Red dotted line
  labs(
    x = "Chromosome",
    y = "Fst"
  ) +
  theme_minimal(base_size = 20) + # Base font size set to 20
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Customize x-axis text
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_blank(),  # Remove plot border
    axis.line = element_line(color = "black"),  # Add main x and y axes
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = chrom_offsets$cumulative_offset + chrom_offsets$offset / 2,
    labels = chrom_offsets$CHROM
  ) +
  scale_y_continuous(
    limits = c(0, max(filtered_df$MeanY, na.rm = TRUE)),
    breaks = c(seq(0, max(filtered_df$MeanY, na.rm = TRUE), by = 0.1), round(top_1_percent_threshold, 2)),
    labels = function(y) ifelse(
      y == round(top_1_percent_threshold, 2), 
      paste0("1% threshold (", y, ")"), 
      y
    )
  )
plot_fst



# Save a copy of the full (unfiltered) dataset
full_df <- filtered_df   # This should be done before you filter for threshold.

# Filter the rows where MeanY exceeds the threshold
filtered_df <- filtered_df[filtered_df$MeanY > top_1_percent_threshold, ]
filtered_df$CHROM <- paste0("chr", filtered_df$CHROM)  # Add "chr" prefix to chromosome numbers

# Order by MeanY (Fst)
ordered_df <- filtered_df[order(-filtered_df$MeanY), ]

# Select rows where MeanY > 0.2
selected_df <- ordered_df[ordered_df$MeanY > top_1_percent_threshold, ]

# Display the result
head(selected_df)


# Load the GFF annotation file
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")

# Filter for genes only
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Load your reference GFF file - replace with the actual path
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Convert Fst data into GRanges
fst_GRanges <- GRanges(
  seqnames = paste0(selected_df$CHROM),
  ranges = IRanges(start = selected_df$WindowStart, end = selected_df$WindowStop),
  Fst = selected_df$MeanY
)

threshold = top_1_percent_threshold

#threshold = 0.2
high_fst_snps <- fst_GRanges[fst_GRanges$Fst >= threshold]

# Find overlapping or nearby genes within 50,000 bp of high Fst SNPs
overlapping_genes <- findOverlaps(high_fst_snps, genes_only)

# Extract gene information
genes_near_high_fst <- genes_only[subjectHits(overlapping_genes)]
gene_names <- mcols(genes_near_high_fst)$ID  # Replace 'ID' with the appropriate column for gene IDs

unique(gene_names)

# Save the results
write.csv(unique(gene_names), "Fst_Spline_5perc_Navarro_maize_NHx_filter_100000_smooth_new.csv", row.names = FALSE)

write.csv(unique(gene_names), "Fst_Spline_5perc_Lasky_sorghum_OlsenP_filter_100000_smooth_new.csv", row.names = FALSE)
getwd()
