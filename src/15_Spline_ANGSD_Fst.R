## ─────────────────────────────────────────────────────────────
## 0) Setup
## ─────────────────────────────────────────────────────────────
# install.packages(c("vroom","dplyr","ggplot2","GenWin","rtracklayer","GenomicRanges","ggrepel","viridis","patchwork"))
library(vroom)
library(dplyr)
library(ggplot2)
library(GenWin)
library(rtracklayer)
library(GenomicRanges)
#library(ggrepel)
#library(viridis)
#library(patchwork)

## >>>> EDIT THESE TWO PATHS <<<<
in_dir  <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Indian_chief/individual"
gff_maize <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

## ─────────────────────────────────────────────────────────────
## 1) Read all per-site Fst (ANGSD) files (gz or plain)
## ─────────────────────────────────────────────────────────────
fns <- list.files(in_dir, pattern = "perSite\\.tsv(\\.gz)?$", full.names = TRUE)
stopifnot(length(fns) > 0)

read_one <- function(fp) {
  df <- vroom::vroom(fp, delim = "\t", show_col_types = FALSE)
  # Normalize column names
  nms <- tolower(names(df))
  if (grepl("^#chr$", nms[1])) names(df)[1] <- "chr"
  names(df) <- tolower(names(df))
  need <- c("chr","pos","alpha","beta","fst")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in ", basename(fp), ": ", paste(miss, collapse=","))
  df %>%
    mutate(chr = as.character(chr)) %>%
    filter(is.finite(fst), fst >= 0, fst <= 1) %>%
    arrange(pos)
}

# read & combine (but we’ll still spline per chr)
fns <- fns[1]
all_list <- lapply(fns, read_one)
all_per_site <- bind_rows(all_list)

# Keep only chr1..chr10; your files are already “chr1” style
all_per_site <- all_per_site %>%
  #filter(chr %in% paste0("chr", 1:10))
  filter(chr %in% paste0("chr", 1))

## ─────────────────────────────────────────────────────────────
## 2) Spline per chromosome (GenWin on per-site Fst)
##    Tune 'smoothness' if you want more/less smoothing
## ─────────────────────────────────────────────────────────────
#chromosomes <- paste0("chr", 1:10)
chromosomes <- paste0("chr", 1)
smooth_kb   <- 100000  # try 100 kb first; you can bump to 200–500 kb
chr_windows <- list()

for (ch in chromosomes) {
  x <- all_per_site %>% filter(chr == ch)
  if (nrow(x) < 100) next  # not enough points to spline
  # ensure strictly increasing map (positions)
  x <- x[order(x$pos), ]
  # splineAnalyze on FST ~ pos
  sp <- GenWin::splineAnalyze(
    Y           = x$fst,
    map         = x$pos,
    smoothness  = smooth_kb,   # in bp
    plotRaw     = FALSE,
    plotWindows = FALSE,
    method      = 4
  )
  
  w <- sp$windowData
  # add chr + keep clean columns
  w$CHROM <- ch
  chr_windows[[ch]] <- w
}

combined_df <- bind_rows(chr_windows)
stopifnot(nrow(combined_df) > 0)

# Clean types & order
combined_df <- combined_df %>%
  mutate(
    CHROM       = factor(CHROM, levels = paste0("chr", 1:10)),
    WindowStart = as.numeric(WindowStart),
    WindowStop  = as.numeric(WindowStop)
  ) %>%
  arrange(CHROM, WindowStart)

# Drop windows with no SNPs or NA
combined_df <- combined_df %>%
  filter(SNPcount >= 1) %>%
  filter(complete.cases(.))

## ─────────────────────────────────────────────────────────────
## 3) Build cumulative chromosome coordinates for plotting
## ─────────────────────────────────────────────────────────────
chrom_offsets <- combined_df %>%
  group_by(CHROM) %>%
  summarise(
    chr_len = max(WindowStop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cumulative_offset = lag(cumsum(chr_len), default = 0),
    chr_mid = cumulative_offset + chr_len / 2
  )

combined_df <- combined_df %>%
  left_join(chrom_offsets, by = "CHROM") %>%
  mutate(
    mid         = (WindowStart + WindowStop) / 2,
    cumulative_pos = mid + cumulative_offset
  )

str(combined_df)
saveRDS(combined_df, file = file.path(in_dir, "FST_spline_windows_combined_chr1.rds"))

## ─────────────────────────────────────────────────────────────
## 4) Threshold (top 1%) and plots
## ─────────────────────────────────────────────────────────────
top1 <- quantile(combined_df$MeanY, 0.99, na.rm = TRUE)

# Manhattan-style plot of smoothed FST
p1 <- ggplot(combined_df, aes(x = cumulative_pos / 1e6, y = MeanY)) +
  geom_line(aes(group = CHROM),
            linewidth = 0.55, alpha = 0.9, colour = "#2E86AB") +
  geom_hline(yintercept = top1, linetype = "dotdash",
             colour = "#E64B35", linewidth = 0.8) +
  scale_x_continuous("Genomic position (Mb)",
                     breaks = chrom_offsets$chr_mid / 1e6,
                     labels = levels(combined_df$CHROM),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expression(italic(F)[ST])) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave(file.path(in_dir, "FST_spline_manhattan.png"), p1, width = 11, height = 4.5, dpi = 300)

## ─────────────────────────────────────────────────────────────
## 5) Top 1% windows → genes (maize)
## ─────────────────────────────────────────────────────────────
ref_GRanges <- rtracklayer::import(gff_maize)
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# FST windows as GRanges
win_gr <- GRanges(
  seqnames = combined_df$CHROM,
  ranges   = IRanges(start = combined_df$WindowStart, end = combined_df$WindowStop),
  Fst      = combined_df$MeanY
)

high_gr <- win_gr[ which(mcols(win_gr)$Fst >= top1) ]
hits    <- findOverlaps(high_gr, genes_only, ignore.strand = TRUE)

top_with_genes <- cbind(
  as.data.frame(high_gr[queryHits(hits)])[, c("seqnames","start","end","Fst")],
  as.data.frame(genes_only[subjectHits(hits)])[, c("seqnames","start","end","ID","Name","biotype")]
)
names(top_with_genes)[1:3] <- c("chr","win_start","win_end")

# write results
out_csv <- file.path(in_dir, "FST_spline_top1pct_windows_annotated.csv")
write.csv(top_with_genes, out_csv, row.names = FALSE)
message("[OK] wrote: ", out_csv)

## Optional: just unique gene IDs near top windows
unique_gene_ids <- unique(top_with_genes$ID)
write.csv(unique_gene_ids, file.path(in_dir, "FST_spline_top1pct_unique_geneIDs.csv"), row.names = FALSE)
