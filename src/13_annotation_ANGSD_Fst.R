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
### Change the chromosome name and the Chr1 from NC....
################################################################################

# cd /Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window_Fst

# cat > map.txt <<'EOF'
# NC_050096.1 chr1
# NC_050097.1 chr2
# NC_050098.1 chr3
# NC_050099.1 chr4
# NC_050100.1 chr5
# NC_050101.1 chr6
# NC_050102.1 chr7
# NC_050103.1 chr8
# NC_050104.1 chr9
# NC_050105.1 chr10
# EOF

# awk '{gsub(/\./,"\\."); printf "s/%s/%s/g\n",$1,$2}' map.txt > map.sed
# 
# # Edit each TSV in place
# for f in *.tsv; do
# sed -i '' -f map.sed "$f"
# done
# 

# while read -r old new; do
# for f in *"$old"*".tsv"; do
# [[ -e "$f" ]] || continue
# mv "$f" "${f//$old/$new}"
# done
# done < map.txt


################################################################################
### LOAD THE REF FILES
################################################################################

ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


################################################################################
### LOAD THE FST FILES
################################################################################

# Read all window FST tables
indir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window_Fst"
fs <- list.files(indir, pattern = "^chr\\d+_J01_J14\\.win100k_step10k\\.tsv$", full.names = TRUE)

read_one <- function(f) {
  # Your files look like: region, chr, midPos, Nsites, Fst
  # Some ANGSD builds label the last column; some don’t.
  dt <- fread(f)
  # Standardize names
  if (ncol(dt) == 5) {
    setnames(dt, c("region","chr","midPos","Nsites","Fst"))
  } else {
    stop(paste("Unexpected columns in", basename(f)))
  }
  dt$file <- basename(f)
  dt
}

allwin <- rbindlist(lapply(fs, read_one))

# 2) QC / parsing
# If you’d rather compute starts/ends from midPos, do it here.
# You used -win 100000 -step 10000, so half-window = 50,000
half_win <- 50000L
allwin <- allwin %>%
  mutate(win_start = pmax(1L, as.integer(midPos - half_win + 1L)),
         win_end   = as.integer(midPos + half_win))

# 3) Choose cutoff: top 1% genome-wide
# optionally filter to autosomes only if needed; your names are chr1..chr10 already
cutoff <- quantile(allwin$Fst, probs = 0.99, na.rm = TRUE)

# Also drop low-coverage windows (tweak threshold to your taste)
min_sites <- 50000L
outliers <- allwin %>%
  filter(!is.na(Fst), Nsites >= min_sites, Fst >= cutoff)

cat(sprintf("99th%% cutoff = %.6f; outlier windows = %d / %d (%.2f%%)\n",
            cutoff, nrow(outliers), nrow(allwin), 100*nrow(outliers)/nrow(allwin)))

# 4) Make GRanges for outlier windows
win_gr <- GRanges(
  seqnames = outliers$chr,
  ranges   = IRanges(start = outliers$win_start, end = outliers$win_end),
  midPos   = outliers$midPos,
  Nsites   = outliers$Nsites,
  Fst      = outliers$Fst,
  file     = outliers$file
)

# Ensure seqlevels match gene object (both use "chr1", etc.)
# If your gene GFF has different naming, harmonize here.
seqlevelsStyle(win_gr) <- seqlevelsStyle(genes_only)

# 5) Overlap windows with genes
hits <- findOverlaps(win_gr, genes_only, ignore.strand = TRUE)

# Pull a tidy table
gene_meta <- mcols(genes_only)
# pick sane columns if present
gene_df <- data.frame(
  gene_idx = subjectHits(hits),
  gene_id  = as.character(gene_meta$ID[subjectHits(hits)]),
  gene_name= as.character(gene_meta$Name[subjectHits(hits)]),
  biotype  = as.character(gene_meta$biotype[subjectHits(hits)]),
  stringsAsFactors = FALSE
)

win_df <- data.frame(
  win_idx  = queryHits(hits),
  chr      = as.character(seqnames(win_gr)[queryHits(hits)]),
  start    = start(win_gr)[queryHits(hits)],
  end      = end(win_gr)[queryHits(hits)],
  midPos   = mcols(win_gr)$midPos[queryHits(hits)],
  Nsites   = mcols(win_gr)$Nsites[queryHits(hits)],
  Fst      = mcols(win_gr)$Fst[queryHits(hits)],
  file     = mcols(win_gr)$file[queryHits(hits)],
  stringsAsFactors = FALSE
)

anno <- cbind(win_df, gene_df) %>%
  arrange(desc(Fst), chr, start, gene_id)

# 6) Write results
outdir <- file.path(indir, "annotated")
dir.create(outdir, showWarnings = FALSE)

f_out1 <- file.path(outdir, "Fst_windows_outliers_annotated.tsv")
f_out2 <- file.path(outdir, "Fst_windows_all_with_cutoff.tsv")

fwrite(anno, f_out1, sep = "\t")
fwrite(allwin %>% mutate(cutoff_99 = cutoff), f_out2, sep = "\t")

cat("Wrote:\n",
    " - ", f_out1, "\n",
    " - ", f_out2, "\n", sep = "")


