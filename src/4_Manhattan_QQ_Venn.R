################################################################################
### LIBRARIES
################################################################################
library(tidyverse)
library(vroom)
library(data.table)
library(purrr)
library(CMplot)
library(qqman)
library(dplyr)
library(GenomicRanges)
library(stringr)


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
  dplyr::rename(Chr = Chromosome, Pos = Position, P.value = p.value) 

farmcpu_raw$log10 <- -log10(farmcpu_raw$P.value)




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



# Load the data
gbs_seeds <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt") %>% dplyr::select(2,1,3,13)
colnames(gbs_seeds) <- c("SNP", "Chromosome", "Position", "N")
gbs_seeds$Chromosome <- paste0("chr",gbs_seeds$Chromosome )



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
CMplot(
  Pmap = combined_df,

  # Plot type and organization
  plot.type = c("m", "c", "q", "d"),
  multraits = FALSE,
  multracks = FALSE,

  # Color settings
  col = c("#4197d8", "#f8c120", "#413496", "#495226",
          "#d60b6f", "#e66519", "#d581b7", "#83d3ad",
          "#7c162c", "#26755d"),
  signal.col = NULL,
  chr.den.col = "black",

  # Point aesthetics
  pch = 19,
  type = "p",
  signal.pch = 19,
  signal.cex = 1.5,
  cex = c(0.5, 1, 1),
  points.alpha = 100L,

  # Threshold settings
  threshold = NULL,
  threshold.col = "red",
  threshold.lty = 2,
  threshold.lwd = 1,
  amplify = TRUE,

  # Chromosome settings
  band = 1,
  chr.border = FALSE,
  chr.labels.angle = 0,
  chr.pos.max = FALSE,
  bin.size = 1e6,

  # Circle plot settings
  cir.chr = TRUE,
  cir.band = 1,
  cir.chr.h = 1.5,
  cir.axis = TRUE,
  cir.axis.col = "black",
  cir.axis.grid = TRUE,
  r = 0.3,
  outward = FALSE,

  # Axis and labels
  ylab = expression(-log[10](italic(p))),
  ylab.pos = 3,
  xticks.pos = 1,
  axis.cex = 1,
  axis.lwd = 1.5,
  lab.cex = 1.5,
  lab.font = 2,

  # Margins and spacing
  mar = c(3, 6, 3, 3),
  mar.between = 0,
  H = 1.5,

  # Output settings
  file = c("jpg", "pdf", "tiff", "png"),
  dpi = 300,
  file.output = TRUE,
  file.name = "",

  # Additional options
  conf.int = TRUE,
  box = FALSE,
  verbose = TRUE
)

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








################################################################################
### INDIVIDUAL PLOTS
################################################################################


# install.packa# install.packages(c("CMplot","dplyr","GenomicRanges"))
library(CMplot)
library(dplyr)
library(GenomicRanges)

## 1) Use the FULL FarmCPU results (all SNPs), not just top hits
# Example read — adjust to your file path:
# farmcpu_all <- readr::read_tsv("FarmCPU.GWAS.Results.txt")

# For this demo I assume you already have a data.frame with these columns:
# SNP, Chr, Pos, P.value
stopifnot(all(c("SNP","Chr","Pos","P.value") %in% names(farmcpu_raw)))

## 2) CMplot format
cm <- farmcpu_raw %>%
  transmute(
    SNP        = as.character(SNP),
    Chromosome = as.integer(Chr),    # must be numeric 1..10
    Position   = as.integer(Pos),
    trait      = as.numeric(P.value) # p-values for CMplot
  ) %>%
  filter(!is.na(Chromosome), !is.na(Position), !is.na(trait)) %>%
  arrange(Chromosome, Position)

# Safety checks (these are usually why points “disappear”)
stopifnot(nrow(cm) > 1000)                    # Manhattan needs LOTS of points
stopifnot(!any(cm$Chromosome <= 0))
stopifnot(!any(duplicated(cm$SNP)))           # duplicated IDs can be dropped by some tools
# If your Chr are like "chr1", convert before here.

## 3) Build highlight set from genes of interest
goi    <- c("Zm00001eb005710", "Zm00001eb399680", "Zm00001eb115450")
pad_bp <- 25000L

goi_gr <- genes_only[mcols(genes_only)$ID %in% goi]
if (length(goi_gr) == 0L) stop("GOI not found in genes_only$ID")
start(goi_gr) <- pmax(1L, start(goi_gr) - pad_bp)
end(goi_gr)   <- end(goi_gr) + pad_bp

snp_gr <- GRanges(seqnames = paste0("chr", cm$Chromosome),
                  ranges   = IRanges(cm$Position, cm$Position))
ov <- findOverlaps(snp_gr, goi_gr, ignore.strand = TRUE)
hl_snps <- cm$SNP[unique(queryHits(ov))]

# Optional: per-gene styles
gene_cols <- c("Zm00001eb005710"="#1b9e77","Zm00001eb399680"="#d95f02","Zm00001eb115450"="#7570b3")
gene_pch  <- c("Zm00001eb005710"=19, "Zm00001eb399680"=17, "Zm00001eb115450"=15)
gene_cex  <- c("Zm00001eb005710"=1.2,"Zm00001eb399680"=1.1,"Zm00001eb115450"=1.1)

hl_goi <- if (length(ov)) mcols(goi_gr)$ID[subjectHits(ov)] else character(0)
hl_col <- unname(gene_cols[hl_goi])
hl_pch <- unname(gene_pch [hl_goi])
hl_cex <- unname(gene_cex [hl_goi])

# Label the lead SNP per gene (min P)
lab <- NULL
if (length(hl_snps)) {
  df_hl <- tibble::tibble(SNP=hl_snps, GOI=hl_goi) %>% left_join(cm[,c("SNP","trait")], by="SNP")
  lab <- df_hl %>% group_by(GOI) %>% slice_min(trait, n=1, with_ties=FALSE) %>% ungroup()
}
hl_text <- if (length(hl_snps)) rep("", length(hl_snps)) else character(0)
if (!is.null(lab) && nrow(lab)) {
  idx <- match(lab$SNP, hl_snps, nomatch = 0)
  hl_text[idx[idx>0]] <- lab$GOI
}

## 4) Whole-genome CMplot
CMplot(cm,
       plot.type      = "m",
       LOG10          = TRUE,                 # tell CMplot these are raw p-values
       #col            = c("#c7dff2","#e6eff7"),
       col            = c("black","lightgrey"),
       threshold      = c(5e-8, 1e-5),
       threshold.col  = c("red","blue"),
       threshold.lty  = c(1,2),
       chr.border     = TRUE,
       amplify        = FALSE,
       highlight      = hl_snps,
       highlight.type = "p",
       highlight.col  = if (length(hl_snps)) hl_col else NULL,
       highlight.pch  = if (length(hl_snps)) hl_pch else NULL,
       highlight.cex  = if (length(hl_snps)) hl_cex else NULL,
       highlight.text = if (length(hl_snps)) hl_text else NULL,
       highlight.text.cex = 0.9,
       file.output    = TRUE,
       file           = "png",
       file.name      = "cmplot_farmcpu_full",
       width          = 14, height = 6, dpi = 300, verbose = TRUE
)

