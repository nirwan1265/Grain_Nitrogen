library(vroom)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# 1) Load & prepare your gene models once
ref_GRanges <- import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# 2) Define chromosomes and trait name
chromosomes <- 1:10
# trait_name  <- "OlsenP_maize"
trait_name  <- "TN_maize"


# 3) Loop over chr1…chr10
for (ch in chromosomes) {
  # a) Read & filter
  fn <- file.path(
    "/Users/nirwantandukar/Downloads/TN_PHG",
    paste0("glm-chr", ch, ".txt")
  )
  dat <- vroom(fn) %>%
    dplyr::select(1, 3, 4, 6) %>%
    filter(p <= 0.05)
  
  # b) Tidy column names
  names(dat) <- c("SNP", "Chromosome", "Position",
                  paste0("pvalue_maize_", trait_name))
  
  # c) Make SNP GRanges, filter & expand
  snp_GR <- GRanges(
    seqnames = Rle(paste0(dat$Chromosome)),
    ranges   = IRanges(dat$Position, width = 1),
    marker   = dat$SNP,
    pvalue   = dat[[4]]
  )
  snp_GR_filt <- snp_GR[mcols(snp_GR)$pvalue <= 0.00000000001]
  snp_GR_exp  <- resize(snp_GR_filt, width = 0, fix = "center")
  
  # d) Find overlaps and build results DF
  hits <- findOverlaps(genes_only, snp_GR_exp)
  odf  <- data.frame(
    Chromosome  = as.character(seqnames(snp_GR_filt[subjectHits(hits)])),
    GeneID       = mcols(genes_only)$ID[ queryHits(hits) ],
    SNP          = mcols(snp_GR_exp)$marker[ subjectHits(hits) ],
    SNP_Position = start(snp_GR_filt[subjectHits(hits)]),
    Gene_Start   = start(genes_only[queryHits(hits)]),
    Gene_End     = end(  genes_only[queryHits(hits)]),
    PValue       = mcols(snp_GR_filt)$pvalue[subjectHits(hits)]
  )
  odf$Relation <- with(odf,
                       ifelse(Gene_Start <= SNP_Position & SNP_Position <= Gene_End,
                              "within",
                              ifelse(SNP_Position < Gene_Start, "upstream", "downstream"))
  )
  
  # e) Write out
  outfn <- paste0("annotation_maize_",
                  trait_name,
                  "_chr", ch,
                  "_GLM.csv")
  write.csv(odf, outfn, row.names = FALSE, quote = FALSE)
}

# 1. List all of the chr1–chr10 output files
files <- list.files(
  path    = ".", 
  pattern = "^annotation_maize_TN_maize_chr[1-9]_GLM\\.csv$|^annotation_maize_TN_maize_chr10_GLM\\.csv$",
  full.names = TRUE
)

# 2. Read and bind them all together
combined_df <- 
  lapply(files, read.csv, stringsAsFactors = FALSE) %>% 
  bind_rows()

unique(combined_df$Chromosome)

# 3. Write out the master file
write.csv(
  combined_df,
  "annotation_maize_TN_allchr_GLM.csv",
  row.names = FALSE,
  quote     = FALSE
)


# Total SNPS
chr1 = 10739387
chr2 = 7244596
chr3 = 8032893
chr4 = 7927187
chr5 = 7673824
chr6 = 5806526
chr7 = 5863869
chr8 = 5587870
chr9 = 5303901
chr10 = 5138391

  
total_snps <-   chr1 + chr2 + chr3 + chr4 + chr5 + chr6 + chr7 + chr8 + chr9 + chr10

# 0.05 significance
alpha <- 0.05

# Bonferroni corrected threshold
bonferroni_threshold <- alpha / total_snps
bonferroni_threshold
#> [1] 7.210598e-10

# In -log10 scale (for Manhattan plots etc.)
log10_thresh <- -log10(bonferroni_threshold)
log10_thresh
#> [1] 9.141699
#> 
#> 
#> 

library(CMplot)
data(cattle50K)
str(dat)

str(cattle50K)
library(stringr)

combined_df_gwas <- combined_df[,c(3,1,4,7)]
colnames(combined_df_gwas) <- c("SNP", "Chromosome", "Position", "pvalue_NHx")
min(combined_df_gwas$pvalue_NHx)

combined_df_gwas <- combined_df_gwas %>%
  mutate(
    SNP = as.factor(SNP),
    chr = as.integer(str_remove(Chromosome, "chr")),  # remove “chr” prefix
    pos = as.integer(Position),
    pval = pvalue_TN
  ) %>%
  dplyr::select(SNP, chr, pos, pval)
str(combined_df_gwas)
-log10(0.00000000000000000001)

CMplot(combined_df_gwas, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(10,70), threshold=c(0.0000000000000000000000001,0.0000000000000000000000001),
       threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(0.5,0.5),signal.pch=c(19,19),
       file="jpg",file.name="TN_maize_gwas",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)









############## OLSEN-P
library(vroom)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# 1) Load & prepare your gene models once
ref_GRanges <- import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# 2) Define chromosomes and trait name
chromosomes <- 1:10
# trait_name  <- "OlsenP_maize"
trait_name  <- "OlsenP_maize"


# 3) Loop over chr1…chr10
for (ch in chromosomes) {
  # a) Read & filter
  fn <- file.path(
    "/Users/nirwantandukar/Downloads/OlsenP_PHG",
    paste0("glm-chr", ch, ".txt")
  )
  dat <- vroom(fn) %>%
    select(1, 3, 4, 6) %>%
    filter(p <= 0.05)
  
  # b) Tidy column names
  names(dat) <- c("SNP", "Chromosome", "Position",
                  paste0("pvalue_maize_", trait_name))
  
  # c) Make SNP GRanges, filter & expand
  snp_GR <- GRanges(
    seqnames = Rle(paste0(dat$Chromosome)),
    ranges   = IRanges(dat$Position, width = 1),
    marker   = dat$SNP,
    pvalue   = dat[[4]]
  )
  snp_GR_filt <- snp_GR[mcols(snp_GR)$pvalue <= 1e-120]
  
  snp_GR_exp  <- resize(snp_GR_filt, width = 5000, fix = "center")
  
  # d) Find overlaps and build results DF
  hits <- findOverlaps(genes_only, snp_GR_exp)
  odf  <- data.frame(
    Chromosome  = as.character(seqnames(snp_GR_filt[subjectHits(hits)])),
    GeneID       = mcols(genes_only)$ID[ queryHits(hits) ],
    SNP          = mcols(snp_GR_exp)$marker[ subjectHits(hits) ],
    SNP_Position = start(snp_GR_filt[subjectHits(hits)]),
    Gene_Start   = start(genes_only[queryHits(hits)]),
    Gene_End     = end(  genes_only[queryHits(hits)]),
    PValue       = mcols(snp_GR_filt)$pvalue[subjectHits(hits)]
  )
  odf$Relation <- with(odf,
                       ifelse(Gene_Start <= SNP_Position & SNP_Position <= Gene_End,
                              "within",
                              ifelse(SNP_Position < Gene_Start, "upstream", "downstream"))
  )
  
  # e) Write out
  outfn <- paste0("annotation_maize_",
                  trait_name,
                  "_chr", ch,
                  "_GLM.csv")
  write.csv(odf, outfn, row.names = FALSE, quote = FALSE)
}
-log10(1e-120)

# 1. List all of the chr1–chr10 output files
files <- list.files(
  path    = ".", 
  pattern = "^annotation_maize_OlsenP_maize_chr[1-9]_GLM\\.csv$|^annotation_maize_OlsenP_maize_chr10_GLM\\.csv$",
  full.names = TRUE
)

# 2. Read and bind them all together
combined_df <- 
  lapply(files, read.csv, stringsAsFactors = FALSE) %>% 
  bind_rows()

# 3. Write out the master file
write.csv(
  combined_df,
  "annotation_maize_OlsenP_allchr_GLM_logp_120.csv",
  row.names = FALSE,
  quote     = FALSE
)
x=1e-20 +1e-22

# Total SNPS
chr1 = 10739387
chr2 = 7244596
chr3 = 8032893
chr4 = 7927187
chr5 = 7673824
chr6 = 5806526
chr7 = 5863869
chr8 = 5587870
chr9 = 5303901
chr10 = 5138391


total_snps <-   chr1 + chr2 + chr3 + chr4 + chr5 + chr6 + chr7 + chr8 + chr9 + chr10

# 0.05 significance
alpha <- 0.05

# Bonferroni corrected threshold
bonferroni_threshold <- alpha / total_snps
bonferroni_threshold
#> [1] 7.210598e-10

# In -log10 scale (for Manhattan plots etc.)
log10_thresh <- -log10(bonferroni_threshold)
log10_thresh
#> [1] 9.141699
#> 
#> 
#> 

library(CMplot)
data(cattle50K)
str(dat)

str(cattle50K)
library(stringr)

combined_df_gwas <- combined_df[,c(3,1,4,7)]
colnames(combined_df_gwas) <- c("SNP", "Chromosome", "Position", "pvalue_OlsenP")

combined_df_gwas <- combined_df_gwas %>%
  mutate(
    SNP = as.factor(SNP),
    chr = as.integer(str_remove(Chromosome, "chr")),  # remove “chr” prefix
    pos = as.integer(Position),
    pval = pvalue_OlsenP
  ) %>%
  dplyr::select(SNP, chr, pos, pval)
str(combined_df_gwas)

 
CMplot(combined_df_gwas, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(2,300), threshold=c(1e-70,1e-70),
       threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(0.5,0.5),signal.pch=c(19,19),
       file="jpg",file.name="TN_maize_gwas",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


