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
gbs_seeds$Chromosome <- paste0("chr",gbs_seeds$Chromosome )


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

