################################################################################
### LIBRARIES
################################################################################
library(tidyverse)
library(vroom)
library(data.table)
library(purrr)
library(CMplot)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)



################################################################################
### GETTING RAW RESULTS FOR MLM, MLMM, BLINK
################################################################################
mlm_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt", 
                 col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, 'H&B.P.Value') %>% dplyr::filter(
                   `H&B.P.Value` <= 0.05
                 )
colnames(mlm_raw)[4] <- "P.value"

BLINK_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_BLINK_3PC_N.txt", 
                   col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value) %>% dplyr::filter(
                     P.value <= 0.0000001
                   )

MLMM_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLMM_3PC_N.txt",
                  col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, 'H&B.P.Value') %>% dplyr::filter(
                    `H&B.P.Value` <= 0.05
                  )
colnames(MLMM_raw)[4] <- "P.value"

## GLM PHG
# 2) Define chromosomes and trait name
chromosomes <- 1:10
# trait_name  <- "OlsenP_maize"
trait_name  <- "TN_maize"


# 3) Loop over chr1…chr10
GLM_raw <- map_dfr(chromosomes, function(ch) {
  fn <- file.path("/Users/nirwantandukar/Downloads/TN_PHG",
                  paste0("glm-chr", ch, ".txt"))
  
  vroom(fn, show_col_types = FALSE) %>%        # fast read
    select(1, 3, 4, 6) %>%                     # keep SNP, Chr, Pos, p
    filter(p <= 1e-25) %>%                     # your p‑value cut‑off
    setNames(c("SNP", "Chr", "Pos", "P.value"))# tidy column names
})

GLM_raw$Chr <- str_remove(GLM_raw$Chr, "chr")

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
  dplyr::rename(Chr = Chromosome, Pos = Position, P.value = p.value) %>%
  dplyr::filter(P.value <= 0.0000001)


colnames(farmcpu_raw)



################################################################################
### LOAD THE REF FILES
################################################################################

ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


################################################################################
### ANNOTATE THE FILES
################################################################################

### USE THE PVALUE HERE FOR DOWNSTREAM PROCESS
# Annotate all in one go
models <- c("MLM", "MLMM", "BLINK", "FarmCPU","GLM")
models <- c("MLM", "MLMM", "BLINK", "FarmCPU")
dfs    <- list(mlm_raw, MLMM_raw, BLINK_raw, farmcpu_raw, GLM_raw)
dfs    <- list(mlm_raw, MLMM_raw, BLINK_raw, farmcpu_raw)
names(dfs) <- models

all_ann <- lapply(models, function(m) {
  df <- dfs[[m]]
  snps <- GRanges(
    seqnames = Rle(paste0("chr", df$Chr)),
    ranges   = IRanges(df$Pos, df$Pos)
  )
  mcols(snps)$marker <- df$SNP
  mcols(snps)$pvalue <- df$P.value
  
  # **NA-safe filtering**
  pv     <- mcols(snps)$pvalue
  #keep   <- !is.na(pv) & pv <= 0.05
  keep   <- !is.na(pv) & pv <= 1
  #keep   <- !is.na(pv) & pv <= 0.0000001
  snp_f  <- snps[keep]
  if (length(snp_f)==0) return(NULL)
  
  snp_e  <- resize(snp_f, width = 25001, fix = "center")
  ol     <- findOverlaps(genes_only, snp_e)
  if (length(ol)==0) return(NULL)
  hg     <- queryHits(ol)
  hs     <- subjectHits(ol)
  
  tibble(
    Model    = m,
    GeneID   = mcols(genes_only)$ID[hg],
    SNP      = mcols(snp_e)$marker[hs],
    Pos      = start(snp_e)[hs],
    P.value  = mcols(snp_e)$pvalue[hs],
    Relation = ifelse(
      start(snp_e)[hs] >= start(genes_only)[hg] &
        start(snp_e)[hs] <= end(genes_only)[hg],
      "within",
      ifelse(
        start(snp_e)[hs] < start(genes_only)[hg],
        "upstream", "downstream"
      )
    )
  )
})

annotated_df <- bind_rows(all_ann)
# annotated_df_farmcpu <- annotated_df %>%
#   dplyr::filter(Model == "FarmCPU")
# write.csv(annotated_df_farmcpu, 
#           "FarmCPU_Annotated_GWAS_Results.csv", 
#           row.names = FALSE)
unique(annotated_df$Model)

quartz()
hist(annotated_df$P.value, breaks = 100, main = "Histogram of P-values", xlab = "P-value", ylab = "Frequency")

table(annotated_df$Model)

library(dplyr)

# 1. Compute max P.value per GeneID × Model
gene_model_max <- annotated_df %>%
  group_by(GeneID, Model) %>%
  summarise(
    max_p = max(P.value, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Define the order you want the models to appear
model_order <- c("MLM", "MLMM", "BLINK","FarmCPU", "GLM")
model_order <- c("MLM", "MLMM", "BLINK","FarmCPU")

# 3. Collapse for each gene
final_summary <- gene_model_max %>%
  # enforce ordering by factor
  mutate(Model = factor(Model, levels = model_order)) %>%
  arrange(GeneID, Model) %>%
  group_by(GeneID) %>%
  summarise(
    Highest_Pvalues = paste(round(max_p, 6), collapse = ";"),
    Model           = paste(as.character(Model), collapse = ", "),
    .groups = "drop"
  )

# 4. Inspect
head(final_summary)


final_summary_sorted <- final_summary %>%
  # 1) Extract the minimum p-value from the semicolon list
  rowwise() %>%
  mutate(
    smallest_p = min(
      as.numeric(strsplit(Highest_Pvalues, ";")[[1]]),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  # 2) Arrange ascending by that smallest p-value
  arrange(smallest_p) %>%
  # 3) Drop the helper column
  select(-smallest_p)

# View the top genes
head(final_summary_sorted)

library(dplyr)

final_summary_sorted <- final_summary_sorted %>%
  mutate(
    No_of_models = lengths(strsplit(Model, ",\\s*"))
  )

# Preview
head(final_summary_sorted)

# Save the final summary to a CSV file
write.csv(final_summary_sorted, 
          "Final_Summary_GWAS_Results_MLM_MLMM_0.05_BLINK_FarmCPU_7.csv", 
          row.names = FALSE)


getwd()




### VENN DIAGRAM

# Define p-value threshold
p_thresh <- 7

# Filter significant genes per model
mlm_genes     <- annotated_df %>% filter(Model == "MLM", P.value < p_thresh)     %>% pull(GeneID) %>% unique()
mlmm_genes    <- annotated_df %>% filter(Model == "MLMM", P.value < p_thresh)    %>% pull(GeneID) %>% unique()
blink_genes   <- annotated_df %>% filter(Model == "BLINK", P.value < p_thresh)   %>% pull(GeneID) %>% unique()
farmcpu_genes <- annotated_df %>% filter(Model == "FarmCPU", P.value < p_thresh) %>% pull(GeneID) %>% unique()
glm_genes     <- annotated_df %>% filter(Model == "GLM", P.value < p_thresh)     %>% pull(GeneID) %>% unique()

# Create a named list for Venn input
venn_input <- list(
  MLM     = mlm_genes,
  MLMM    = mlmm_genes,
  BLINK   = blink_genes,
  FarmCPU = farmcpu_genes,
  GLM     = glm_genes
)

# Generate Venn diagram object
library(VennDiagram)
venn.plot <- venn.diagram(
  x = venn_input,
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "goldenrod", "purple", "darkgreen", "red"),
  alpha = 0.65,
  lwd = 2,
  lty = "solid",
  col = "white",
  cex = 1.4,
  fontface = "bold",
  cat.cex = 1.3,
  cat.fontface = "bold",
  margin = 0.08,
  height = 6,
  width = 6,
  rotation.degree = 0,
  main = "Significant Genes at p < 1e-7",
  main.cex = 1.5,
  main.pos = c(0.5, 1.05)
)

# Plot the Venn diagram
quartz()
venn.plot
grid.newpage()
grid.draw(venn.plot)

# Optional: save as PNG
png("gene_venn_diagram_logp_7.png", width = 5, height = 5, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
