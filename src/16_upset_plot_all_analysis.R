# install.packages(c("readr","UpSetR"))  # if needed
library(readr)
library(UpSetR)


df <- read_csv("results/all_candidate_genes.csv", show_col_types = FALSE)

# compute only combos with size > 0
set_names <- names(sets)
comb_list <- unlist(lapply(seq_along(set_names), function(k)
  combn(set_names, k, simplify = FALSE)), recursive = FALSE)
non_empty <- Filter(function(v) length(Reduce(intersect, sets[v])) > 0, comb_list)

up_df2 <- fromList(sets)
up_df2 <- up_df2[rowSums(up_df2) > 0, , drop = FALSE]

png("UpSet_genes_all_analyses.png", 1600, 900, res = 200)
upset(up_df2,
      intersections       = non_empty,     # <- forces only these combos
      order.by            = "freq",
      sets.bar.color      = "#6baed6",
      main.bar.color      = "grey20",
      matrix.color        = "grey20")
dev.off()



### Get the combinations
# 1) read
df <- read_csv("results/all_candidate_genes.csv", show_col_types = FALSE)

# 2) define your set columns (edit if names differ)
set_cols <- c("Fst_landraces_N",
              "Fst_IC_stover",
              "Fst_J_grain",
              "GWAS_PHG_landraces_N",
              "GWAS_GBS_landraces_N")

# 3) build sets: unique gene IDs per column
sets <- lapply(df[set_cols], function(x){
  x <- trimws(as.character(x))
  unique(na.omit(x))
})
names(sets) <- set_cols

# quick sanity
#sapply(sets, length)

# 4) all non-empty combinations and their intersections
all_combos <- unlist(lapply(seq_along(sets), function(k)
  combn(names(sets), k, simplify = FALSE)), recursive = FALSE)

intersections <- map_df(all_combos, function(nms){
  genes <- Reduce(intersect, sets[nms])
  tibble(
    combo = paste(nms, collapse = " & "),
    k     = length(nms),
    size  = length(genes),
    genes = list(sort(genes))
  )
}) %>%
  filter(size > 0) %>%
  arrange(desc(k), desc(size), combo)

# 5) save a summary CSV (genes column is semicolon-joined for readability)
summary_out <- intersections %>%
  mutate(genes = vapply(genes, function(g) paste(g, collapse = ";"), character(1)))

write.csv(summary_out, "gene_intersections_summary.csv", row.names = F)
