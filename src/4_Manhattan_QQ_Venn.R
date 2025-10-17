################################################################################
### LIBRARIES
################################################################################
library(tidyverse)
library(vroom)
library(data.table)
library(purrr)
library(CMplot)

################################################################################
### GETTING RAW RESULTS FOR MLM, MLMM, BLINK
################################################################################
mlm_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt", 
                 col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

BLINK_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_BLINK_3PC_N.txt", 
                 col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

MLMM_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLMM_3PC_N.txt",
                 col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

################################################################################
### COMBINING THE FARMCPU RESULTS
################################################################################

# Initialize empty list to hold per-chromosome results
gwas_list <- list()

# Loop through chr1 to chr10
for (chr in 1:10) {
  # Construct the file path
  file_path <- paste0("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/FarmCPU/FarmCPU_TN_3PC_maize_chr", chr, ".rds")
  
  # Read the RDS file
  x <- readRDS(file_path)
  
  # Extract the GWAS results table for TN_maize
  gwas_chr <- x$TN_maize$GWAS
  
  # Append to the list
  gwas_list[[chr]] <- gwas_chr
}

# Combine all into a single data frame
farmcpu_raw <- do.call(rbind, gwas_list)

# Optional: check result
farmcpu_raw <- farmcpu_raw %>%
  dplyr::select(SNP, Chromosome, Position, p.value) %>%
  #rename the columns to match the others
  rename(Chr = Chromosome, Pos = Position, P.value = p.value) 



################################################################################
### COMBINING THE RESULTS
################################################################################

# First, rename columns to standard format for merging
mlm_df     <- mlm_raw     |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, MLM     = P.value)
mlmm_df    <- MLMM_raw    |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, MLMM    = P.value)
blink_df   <- BLINK_raw   |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, BLINK   = P.value)
farmcpu_df <- farmcpu_raw |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, FarmCPU = P.value)

# Merge all by SNP, Chromosome, and Position
combined_df <- mlm_df %>%
  full_join(mlmm_df,   by = c("SNP", "Chromosome", "Position")) %>%
  full_join(blink_df,  by = c("SNP", "Chromosome", "Position")) %>%
  full_join(farmcpu_df, by = c("SNP", "Chromosome", "Position"))

# Reorder columns: SNP, Chromosome, Position, then models
combined_df <- combined_df |> dplyr::select(SNP, Chromosome, Position, MLM, MLMM, BLINK, FarmCPU)

?CMplot()

# PLOT MANHATTAN
CMplot(
  Pmap = combined_df,
  plot.type = "m",
  multraits = TRUE,
  
  # Color settings
  col = "grey",
  signal.col = c("dodgerblue", "goldenrod", "purple", "olivedrab3"),
  #signal.col = viridis::viridis(3, option = "D", end = 0.9),  # "D" = viridis, "C" = magma
  threshold.col = c("black", "grey"),
  
  # Point aesthetics
  signal.pch = c(19, 6, 4, 10),
  signal.cex = 0.8,
  points.alpha = 255,
  cex = 0.8,
  
  # Threshold settings
  threshold = c(1e-7,1e-5),
  threshold.lty = 1,
  threshold.lwd = c(1, 1),
  amplify = TRUE,
  
  # Chromosome density
  chr.den.col = NULL,
  axis.lwd = 1,
  
  # Output settings
  file = "jpg",
  file.name = NULL,
  dpi = 300,
  file.output = TRUE,
  width = 12,
  height = 5,
  
  # Legend
  legend.ncol = 3,
  legend.pos = "middle",
  
  # Additional
  verbose = TRUE
)

#### Individual plots
# CMplot(
#   Pmap = combined_df,
#   
#   # Plot type and organization
#   plot.type = c("m", "c", "q", "d"),
#   multraits = FALSE,
#   multracks = FALSE,
#   
#   # Color settings
#   col = c("#4197d8", "#f8c120", "#413496", "#495226", 
#           "#d60b6f", "#e66519", "#d581b7", "#83d3ad", 
#           "#7c162c", "#26755d"),
#   signal.col = NULL,
#   chr.den.col = "black",
#   
#   # Point aesthetics
#   pch = 19,
#   type = "p",
#   signal.pch = 19,
#   signal.cex = 1.5,
#   cex = c(0.5, 1, 1),
#   points.alpha = 100L,
#   
#   # Threshold settings
#   threshold = NULL,
#   threshold.col = "red",
#   threshold.lty = 2,
#   threshold.lwd = 1,
#   amplify = TRUE,
#   
#   # Chromosome settings
#   band = 1,
#   chr.border = FALSE,
#   chr.labels.angle = 0,
#   chr.pos.max = FALSE,
#   bin.size = 1e6,
#   
#   # Circle plot settings
#   cir.chr = TRUE,
#   cir.band = 1,
#   cir.chr.h = 1.5,
#   cir.axis = TRUE,
#   cir.axis.col = "black",
#   cir.axis.grid = TRUE,
#   r = 0.3,
#   outward = FALSE,
#   
#   # Axis and labels
#   ylab = expression(-log[10](italic(p))),
#   ylab.pos = 3,
#   xticks.pos = 1,
#   axis.cex = 1,
#   axis.lwd = 1.5,
#   lab.cex = 1.5,
#   lab.font = 2,
#   
#   # Margins and spacing
#   mar = c(3, 6, 3, 3),
#   mar.between = 0,
#   H = 1.5,
#   
#   # Output settings
#   file = c("jpg", "pdf", "tiff", "png"),
#   dpi = 300,
#   file.output = TRUE,
#   file.name = "",
#   
#   # Additional options
#   conf.int = TRUE,
#   box = FALSE,
#   verbose = TRUE
# )

### QQPLOT

CMplot(
  Pmap = combined_df,
  plot.type = "q",
  multraits = TRUE,
  
  # Color settings
  col = c("dodgerblue", "goldenrod", "purple", "olivedrab3"),
  signal.col = "red",
  conf.int.col = NULL,
  
  # Point aesthetics
  signal.pch = c(19, 6, 4, 10),
  signal.cex = 0.8,
  points.alpha = 255,
  cex = 0.8,
  
  # Threshold settings
  threshold = 1e-7,
  threshold.col = "red",
  threshold.lty = 2,
  threshold.lwd = 1,
  
  # Axis and labels
  ylab = expression(-log[10](italic(p))),
  ylab.pos = 2,
  ylim = c(0, 16),
  axis.cex = 1,
  axis.lwd = 1,
  
  # Confidence intervals
  conf.int = TRUE,
  
  # Output settings
  file = "jpg",
  dpi = 300,
  #width = 5,
  #height = 5,
  file.output = TRUE,
  file.name = NULL,
  
  # Layout
  box = FALSE,
  verbose = TRUE
)




#### VENN DIAGRAM
# Load required package
library(VennDiagram)

# Define threshold
p_thresh <- 1e-7

# Filter significant SNPs
mlm_sig     <- na.omit(mlm_raw$SNP     [mlm_raw$P.value     < p_thresh])
mlmm_sig    <- na.omit(MLMM_raw$SNP    [MLMM_raw$P.value    < p_thresh])
blink_sig   <- na.omit(BLINK_raw$SNP   [BLINK_raw$P.value   < p_thresh])
farmcpu_sig <- na.omit(farmcpu_raw$SNP[farmcpu_raw$P.value < p_thresh])


# Create named list for input
venn_input <- list(
  MLM     = mlm_sig,
  MLMM    = mlmm_sig,
  BLINK   = blink_sig,
  FarmCPU = farmcpu_sig
)



venn.plot <- venn.diagram(
  x = venn_input,
  filename = NULL,  # Draw to R object instead of file
  output = TRUE,    # Return grob object
  
  # Color settings
  fill = c("dodgerblue", "goldenrod", "purple", "darkgreen"),  # Matplotlib tableau colors
  alpha = 0.65,     # Slightly more transparency for better overlap visibility
  
  # Circle borders
  lwd = 2,          # Thicker border lines
  lty = "solid",    # Solid line style
  col = "white",    # White borders for clean look
  
  # Numeric labels
  cex = 1.4,        # Larger size for numbers
  fontface = "bold",
  #fontfamily = "sans",  # Clean font
  
  # Category labels
  cat.cex = 1.3,
  cat.fontface = "bold",
  #cat.fontfamily = "sans",
  #cat.pos = c(-30, 30, 180, 0),  # Better angular positions
  #cat.dist = c(0.1, 0.1, 0.1, 0.1),  # Distance from circles
  
  # Margins and scaling
  margin = 0.08,
  height = 6,       # Inches
  width = 6,
  
  # Rotation (if needed)
  rotation.degree = 0,
  
  # Special effects
  # ext.text = TRUE,  # Better text positioning
  # ext.line.lwd = 1,
  # ext.dist = -0.15,
  # ext.length = 0.8,
  
  # Title (optional)
  main = "SNPs at -log10(p) = 7",
  main.cex = 1.5,
  main.pos = c(0.5, 1.05)
)

# To display the plot
quartz()
grid.draw(venn.plot)

# To save as high-quality image
# Save as PNG
png("venn_diagram.png", width=5, height=5, units="in", res=300)
grid.draw(venn.plot)
dev.off()


?GAPIT()
