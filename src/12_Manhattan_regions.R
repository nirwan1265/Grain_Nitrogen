library(CMplot)
library(vroom)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
# data(pig60K)
# head(pig60K)

# Load the data
gbs_seeds <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt") %>% dplyr::select(2,1,3,13)
colnames(gbs_seeds) <- c("SNP", "Chromosome", "Position", "N")
gbs_seeds$Chromosome <- paste0(gbs_seeds$Chromosome, "GBS")


# 2) Define chromosomes and trait name
chromosomes <- 1:10
# trait_name  <- "OlsenP_maize"
trait_name  <- "TN_maize"

phg_seeds <- map_dfr(chromosomes, function(ch) {
  fn <- file.path("/Users/nirwantandukar/Downloads/TN_PHG",
                  paste0("glm-chr", ch, ".txt"))
  
  vroom(fn, show_col_types = FALSE) %>%        # fast read
    dplyr::select(1, 3, 4, 6) %>%                     # keep SNP, Chr, Pos, p
    setNames(c("SNP", "Chromosome", "Position", "N"))# tidy column names
})

colnames(phg_seeds) <- c("SNP", "CHR",  "BP","P")
phg_seeds <- phg_seeds[!is.na(phg_seeds$P), ]

# save rds
saveRDS(phg_seeds,"phg_seeds.RDS")
phg_seeds <- readRDS("phg_seeds.RDS")

# Get the snps
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# Define buffer size
buffer <- 1  # 25 kb

# Step 1: Extract gene info
target_gene <- "Zm00001eb421970"
gene_row <- genes_only[mcols(genes_only)$ID == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- phg_seeds[phg_seeds$CHR == paste0("chr",gene_chr), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(BP >= (gene_start - buffer) & BP <= (gene_end + buffer)) %>%
  arrange(BP)

# Final output
snps_in_window

# Remove the string "chr" from CHR column in phg_seeds
phg_seeds$CHR <- as.numeric(gsub("chr", "", phg_seeds$CHR))

# Check the hits:
region_hits <- subset(
  phg_seeds,
  CHR == 10 &
    BP  >= gene_start &
    BP  <= gene_end
)
nrow(region_hits)  

quartz()
manhattan(subset(phg_seeds, CHR == 9), highlight = snps_in_window$SNP, xlim = c(gene_start,gene_end), main = "NMX")

