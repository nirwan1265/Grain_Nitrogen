# install.packages(c("vroom","dplyr","readr","ggplot2","IRanges","GenomicRanges"), Ncp=1)
library(vroom)
library(dplyr)
library(readr)
library(ggplot2)
library(GenomicRanges)

# 1) Read all your window files
indir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window"
files <- list.files(indir, pattern="^chr[0-9]+_IC01_IC14\\.win100k_step10k\\.clean\\.tsv$", full.names = TRUE)

win <- vroom::vroom(files, col_types = "ciiiid") %>%  # chr,start,end,midPos,Nsites,Fst
  mutate(chr = factor(chr, levels=paste0("chr", 1:10))) %>% 
  arrange(chr, start)

# 2) Choose a cutoff (top 1% genome-wide by default)
cutoff <- quantile(win$Fst, probs = 0.99, na.rm = TRUE)
cutoff

win <- win %>% mutate(outlier = Fst >= cutoff)

# 3) Quick genome-wide plot (Manhattan-like using cumulative positions)
#    Build cumulative offsets so chromosomes plot in order
chr_lengths <- win %>% group_by(chr) %>% summarise(chr_len = max(end, na.rm=TRUE))
chr_lengths <- chr_lengths %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))
win <- win %>% left_join(chr_lengths, by="chr") %>%
  mutate(pos_cum = midPos + offset)

# x-axis tick positions (mid of each chromosome)
chr_ticks <- chr_lengths %>%
  mutate(mid = offset + chr_len/2)

plot <- ggplot(win, aes(pos_cum, Fst, color = chr)) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_hline(yintercept = cutoff, linetype = 2) +
  scale_x_continuous(breaks = chr_ticks$mid, labels = as.character(chr_ticks$chr)) +
  guides(color = "none") +
  labs(x = "Chromosome", y = "Fst (100kb windows, 10kb step)",
       title = "Genome-wide sliding-window Fst (IC01 vs IC14)") +
  theme_minimal(base_size = 12)


quartz()
plot
# 4) Convert top windows to GRanges for gene overlap
top_win_gr <- with(win %>% filter(outlier), 
                   GRanges(seqnames = chr, ranges = IRanges(start = start, end = end)))

# References
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

## hits: overlaps between top windows and genes
hits <- findOverlaps(top_win_gr, genes_only, ignore.strand = TRUE)

## Pull fields directly from the query (windows) and subject (genes)
qq <- queryHits(hits)
ss <- subjectHits(hits)

top_with_genes <- tibble(
  chr       = as.character(seqnames(top_win_gr))[qq],
  win_start = start(top_win_gr)[qq],
  win_end   = end(top_win_gr)[qq],
  gene_chr  = as.character(seqnames(genes_only))[ss],
  gene_start= start(genes_only)[ss],
  gene_end  = end(genes_only)[ss],
  gene_id   = mcols(genes_only)$ID[ss],
  gene_name = mcols(genes_only)$Name[ss],
  biotype   = mcols(genes_only)$biotype[ss]
)

## Make sure join keys have the same type as in top_tbl
top_tbl <- top_tbl %>% mutate(chr = as.character(chr))

## Add Fst/Nsites/midPos back onto the overlap rows
top_annot <- top_with_genes %>%
  inner_join(top_tbl, by = c("chr" = "chr",
                             "win_start" = "start",
                             "win_end"   = "end"))

## Write out if you want
write.csv(top_annot,"Fst_top1pct_windows_with_genes.csv", row.names = F, quote=F)
getwd()
## Quick sanity check
dplyr::glimpse(top_annot)